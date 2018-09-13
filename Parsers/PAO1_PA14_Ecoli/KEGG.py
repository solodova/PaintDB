from Bio.KEGG.REST import kegg_list, kegg_get, kegg_conv
from Bio.KEGG.KGML.KGML_parser import read
from Schema import Interactor, Metabolite, OrthologEcoli, Interaction, InteractionReference, InteractionSource

kegg_compounds = {}

# function to parse kegg interactions for PAO1 and PA14
def parse_pseudomonas(session):
    # get all the kegg compounds if they were not obtained yet
    if not kegg_compounds:
        get_kegg_compounds()
    # add PAO1 and PA14 KEGG sources
    # Note: is_experimental is set to 2 since no way to deteremine if detection method was experimental or not
    source_PAO1 = InteractionSource(data_source='KEGG(PAO1)', is_experimental=2)
    source_PA14 = InteractionSource(data_source='KEGG(PA14)', is_experimental=2)
    session.add(source_PAO1), session.add(source_PA14), session.commit()
    # parse PAO1 and PA14 Kegg interactions
    parse_kegg('pae', 'PAO1', 'KEGG(PAO1)', session)
    parse_kegg('pau', 'PA14', 'KEGG(PA14)', session)

# function to parse kegg interactions from Ecoli
def parse_ecoli(session):
    # get all kegg compounds if they were not obtained yet
    if not kegg_compounds:
        get_kegg_compounds()
    # update metabolite info for existing metabolites which may be missing ids
    update_metabolite_info_kegg(session)
    # add Ecoli KEGG source
    # Note: is_experimental is set to 2 since no way to determine if detection method was experimental or not
    source = InteractionSource(data_source='KEGG(Ecoli)', is_experimental=2)
    session.add(source), session.commit()
    # parse kegg interactions from Ecoli, doing PAO1 and PA14 ortholog mapping separately
    parse_kegg('eco', 'PAO1', 'KEGG(Ecoli)', session)
    parse_kegg('eco', 'PA14', 'KEGG(Ecoli)', session)

#function to find type of interaction based on attribute abbreviation given
def find_type_kegg(attrib):
    return {
        'ECrel': 'enzyme-enzyme relation',
        'PPrel': 'protein-protein interaction',
        'GErel': 'gene expression interaction',
        'PCrel': 'protein-compound interaction',
    }[attrib]

# make sure to run this before calling parse_?_KEGG()
def get_kegg_compounds():
    # fill in compounds dictionary (kegg_id: {name: '', pubchem: '', chebi: ''})
    for compound in kegg_list(database='compound'):
        kegg_compounds[compound[4:10]] = {}
        kegg_compounds[compound[4:10]]['name'] = compound.split('\t')[1].split(';')[0].rstrip()
        kegg_compounds[compound[4:10]]['chebi'] = None
        kegg_compounds[compound[4:10]]['pubchem'] = None

    # had to change kegg_conv source code to accept 'chebi' as target db for kegg_conv
    for cpd_id in kegg_conv('chebi', 'compound').read().split('cpd:'):
        if cpd_id != '':
            kegg_compounds[cpd_id[:6]]['chebi'] = cpd_id.split('chebi:')[1].rstrip()
    for cpd_id in kegg_conv('pubchem', 'compound').read().split('cpd:'):
        if cpd_id != '':
            kegg_compounds[cpd_id[:6]]['pubchem'] = cpd_id.split('pubchem:')[1].rstrip()

# function to update metabolite info for existing metabolites with missing ids
def update_metabolite_info_kegg(session):
    for kegg_id in kegg_compounds:
        metabolite = None
        pubchem, name, chebi =  kegg_compounds[kegg_id]['pubchem'],  kegg_compounds[kegg_id]['name'], \
                                kegg_compounds[kegg_id]['chebi']
        if pubchem is not None:
            metabolite = session.query(Metabolite).filter_by(pubchem = pubchem).first()
        if (metabolite is None) & (chebi is not None):
            metabolite = session.query(Metabolite).filter_by(chebi=chebi).first()
        if metabolite is None:
            metabolite = session.query(Metabolite).filter_by(kegg=kegg_id).first()
        if (metabolite is None) & (name is not None):
            metabolite = session.query(Metabolite).filter_by(name=name).first()

        if metabolite is not None:
            if metabolite.kegg is None:
                metabolite.kegg = kegg_id
            if metabolite.pubchem is None:
                metabolite.pubchem = pubchem
            if metabolite.chebi is None:
                metabolite.chebi = chebi
            if metabolite.name is None:
                metabolite.name = name

def parse_kegg(org_id, strain, sourcedb, session):
    # get pathways for organism specified by org_id
    pathways = kegg_list(database='pathway', org=org_id).read().split('path:')
    path_ids = []

    # make list of path ids to iterate through
    for path in pathways:
        if path != '':
            path_ids.append(path[:8])

    # iterate through each path and obtain interactions
    for path in path_ids:
        # get kgml representation of path
        kgml_path = read(kegg_get(path, option='kgml'))
        path_name = kgml_path._getname()
        # dictionary of compounds in current path (node_id: kegg_id)
        #   compound._getid() returns node id (only relevant in context of current path)
        #   compound._getname() returns kegg id (relevant in overall KEGG DB)
        compound_ids = {}
        for compound in kgml_path.compounds:
            compound_ids[compound._getid()] = compound._getname()[-6:]
        # go through each relation in path
        for relation in kgml_path.relations:
            relation_type = relation.element.attrib['type']

            # ignore maplink relations
            if relation_type == 'maplink': continue
            # relation._getentry1/2() returns  protein id (locus) or compound id (KEGG id)
            entries = [relation._getentry1()._getname(), relation._getentry2()._getname()]
            # if one or both interactors are listed as undefined, move on to next interaction
            if (entries[0] == 'undefined') | (entries[1] == 'undefined'): continue
            # list to hold existing interactors
            interactors = [[], []]
            # list to hold new metabolite ids for interactions with metabolites not yet in the database
            new_metabolites = [[], []]
            # go through each entry in the relation
            for num in range(0, 2):
                # each entry may contain >1 id; go through all of them
                for id in entries[num].split(' '):
                    if id == '': continue
                    # if interactor is not protein or compound, continue
                    if (id.split(':')[0] != org_id) & (id.split(':')[1] not in kegg_compounds): continue

                    # check if the id is a kegg id by searching in kegg_compounds
                    kegg_id= None
                    if id.split(':')[1] in kegg_compounds:
                        kegg_id = id.split(':')[1]

                    # check if interactor (protein) already exists
                    if (kegg_id is None) & (org_id != 'eco'):
                        interactor = session.query(Interactor).get(id.split(':')[1])
                        if interactor is not None:
                            # make sure to add None value; this will be needed to create interaction reference later
                            # None is appended rather than the interactor id because the interactor is not an ortholog
                            interactors[num].append([interactor, None])
                    # if it doesnt exist, it's not a valid protein, so check if it is a valid compound
                    elif kegg_id is not None:
                        interactor = session.query(Metabolite).filter_by(kegg = kegg_id).first()
                        # if metabolite with id was not found, append the kegg_id to new_metabolites to create
                        if interactor is None:
                            new_metabolites[num].append(kegg_id)
                        else:
                            # if the metabolite was found, add it to the existing interactor list
                            interactors[num].append([interactor, interactor.id])
                    # if parsing E. coli path, add all orthologs to interactor list
                    elif org_id == 'eco':
                        for ortholog in session.query(OrthologEcoli).filter_by(ortholog_id = id.split(':')[1],
                                                                                strain_protein = strain).all():
                            if ortholog is not None:
                                # add the id of the ecoli protein for the interaction reference later
                                interactors[num].append([ortholog.protein, id.split(':')[1]])

            # create list of interactor pairs from two separate lists
            interactor_pairs = []
            # create interactor pairs from interactors which already exist in db
            for interactor1 in interactors[0]:
                for interactor2 in interactors[1]:
                    if (interactor1[0].type != 'm') | (interactor2[0].type != 'm'):
                        interactor_pairs.append([interactor1, interactor2])
            # create interactor pair from interactors and new metabolites
            for interactor1 in interactors[0]:
                for id in new_metabolites[1]:
                    # ignore interactor pairs which would result in m-m interactions
                    if interactor1[0].type == 'm': continue
                    # Note: can query metabolite with kegg only because we updated the metabolite info first
                    metabolite = session.query(Metabolite).filter_by(kegg = id).first()
                    if metabolite is None:
                        metabolite = Metabolite(id = id, kegg = id, pubchem = kegg_compounds[id]['pubchem'],
                                                chebi = kegg_compounds[id]['chebi'])
                        session.add(metabolite)
                    interactor_pairs.append([interactor1, [metabolite, metabolite.id]])
            for interactor1 in interactors[1]:
                for id in new_metabolites[0]:
                    if interactor1[0].type == 'm': continue
                    metabolite = session.query(Metabolite).filter_by(kegg = id).first()
                    if metabolite is None:
                        metabolite = Metabolite(id = id, kegg = id, pubchem = kegg_compounds[id]['pubchem'],
                                                chebi = kegg_compounds[id]['chebi'])
                        session.add(metabolite)
                    interactor_pairs.append([interactor1, [metabolite, metabolite.id]])

            # if no interactor pairs were found, move on the the next interaction
            if len(interactor_pairs) == 0: continue

            # get all intermediates in reaction of type compound
            intermeds = []
            for subtype in relation.element.iter(tag='subtype'):
                # if the subtype element is a compound, get its node id
                if 'compound' in subtype.attrib:
                    compound_node_id = subtype.attrib['compound']
                    if compound_node_id is None: continue
                    # if the node id was not stored in the compound ids for this path, move on to the next sybtype
                    if int(compound_node_id) not in compound_ids: continue
                    # if compound id is valid, either add existing matching metabolite or create new one and add
                    kegg_id = compound_ids[int(compound_node_id)]
                    metabolite = session.query(Metabolite).filter_by(kegg = kegg_id).first()
                    if metabolite is None:
                        metabolite = Metabolite(id=kegg_id, name=kegg_compounds[kegg_id]['name'],
                                                pubchem=kegg_compounds[kegg_id]['pubchem'],
                                                chebi=kegg_compounds[kegg_id]['chebi'], kegg=kegg_id)
                        session.add(metabolite)
                    intermeds.append([metabolite, metabolite.id])

            # add protein - intermediate interactor pairs
            for interactor_list in interactors:
                for interactor in interactor_list:
                    if interactor[0].type != 'm':
                        for intermed in intermeds:
                            interactor_pairs.append([interactor, intermed])

            # go through each interaction pair and add interaction if it doesnt exist yet
            for interactor_pair in interactor_pairs:
                homogenous = (interactor_pair[0][0] == interactor_pair[1][0])
                interaction = session.query(Interaction).filter(Interaction.interactors.contains(interactor_pair[0][0]),
                                                                Interaction.interactors.contains(interactor_pair[1][0]),
                                                                Interaction.homogenous == homogenous).first()

                source = session.query(InteractionSource).filter_by(data_source=sourcedb).first()
                #create interaction if it doesnt exist yet, add source to its sources if it isn't already
                if interaction is None:
                    interaction = Interaction(type=interactor_pair[0][0].type + '-' + interactor_pair[1][0].type,
                                              strain=strain, homogenous=homogenous,
                                              interactors=[interactor_pair[0][0], interactor_pair[1][0]])
                    interaction.sources.append(source)
                    if org_id == 'eco':
                        interaction.ortholog_derived = 'Ecoli'
                    session.add(interaction), session.commit()
                elif source not in interaction.sources:
                    interaction.sources.append(source)

                # in case the interaction already existed, make sure interactor_a and interactor_b variables for
                # new interaction reference match up with the first and second interactors of the existing
                # interaction (so it's easy to see which Pseudomonas interactor matches up with which Ecoli
                # ortholog if the org id is eco)
                interactor_a, interactor_b = None, None
                if org_id == 'eco':
                    if interaction.interactors[0] == interactor_pair[0][0]:
                        interactor_a = interactor_pair[0][1]
                        interactor_b = interactor_pair[1][1]
                    else:
                        interactor_b = interactor_pair[0][1]
                        interactor_a = interactor_pair[1][1]

                # search for reference
                reference = session.query(InteractionReference).filter_by(source_db='kegg',
                                                                          comment='in ' + path_name + ' path',
                                                                          interactor_a=interactor_a,
                                                                          interactor_b=interactor_b).first()
                # if the reference doesnt exist, create it, add it to the interaction's references and add the source
                # to the reference's sources
                if reference is None:
                    reference = InteractionReference(source_db='kegg', comment='in ' + path_name + ' path',
                                                     interactor_a=interactor_a, interactor_b=interactor_b)
                    interaction.references.append(reference)
                    reference.sources.append(source)
                # if the reference does exist, add it to the interaction's reference list and add the source to the
                # reference's source list if it isn't there already
                else:
                    if interaction not in reference.interactions:
                        reference.interactions.append(interaction)
                    if source not in reference.sources:
                        reference.sources.append(source)

    session.commit()
    print(sourcedb, session.query(Interaction).count())