from Bio.KEGG.REST import kegg_list, kegg_get, kegg_conv
from Bio.KEGG.KGML.KGML_parser import read
from Schema import Interactor, Metabolite, OrthologEcoli, Interaction, InteractionReference, InteractionSource

kegg_compounds = {}

def parse_pseudomonas(session):
    if not kegg_compounds:
        get_kegg_compounds()
    source_PAO1 = InteractionSource(data_source='KEGG(PAO1)', is_experimental=2)
    source_PA14 = InteractionSource(data_source='KEGG(PA14)', is_experimental=2)
    session.add(source_PAO1), session.add(source_PA14), session.commit()
    parse_kegg('pae', 'PAO1', 'KEGG(PAO1)', session)
    parse_kegg('pau', 'PA14', 'KEGG(PA14)', session)

def parse_ecoli(session):
    if not kegg_compounds:
        get_kegg_compounds()
    #update_metabolite_info(session)
    source = InteractionSource(data_source='KEGG(Ecoli)', is_experimental=2)
    session.add(source), session.commit()
    parse_kegg('eco', 'PAO1', 'KEGG(Ecoli)', session)
    parse_kegg('eco', 'PA14', 'KEGG(Ecoli)', session)

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

    # had to change source code to accept 'chebi' as target db for kegg_conv
    for cpd_id in kegg_conv('chebi', 'compound').read().split('cpd:'):
        if cpd_id != '':
            kegg_compounds[cpd_id[:6]]['chebi'] = cpd_id.split('chebi:')[1].rstrip()
    for cpd_id in kegg_conv('pubchem', 'compound').read().split('cpd:'):
        if cpd_id != '':
            kegg_compounds[cpd_id[:6]]['pubchem'] = cpd_id.split('pubchem:')[1].rstrip()


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

    for path in pathways:
        if path != '':
            path_ids.append(path[:8])

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

            if relation_type == 'maplink': continue
            # relation._getentry1/2() returns  protein id (locus) or compound id (KEGG id)
            entries = [relation._getentry1()._getname(), relation._getentry2()._getname()]
            if (entries[0] == 'undefined') | (entries[1] == 'undefined'): continue
            interactors = [[], []]
            new_metabolites = [[], []]
            # go through each entry in the relation, find/create interactors_sif
            for num in range(0, 2):
                # each entry may contain >1 id; go through all of them
                for id in entries[num].split(' '):
                    if id == '': continue
                    # if interactor is not protein or compound, continue
                    if (id.split(':')[0] != org_id) & (id.split(':')[1] not in kegg_compounds): continue

                    kegg_id= None
                    if id.split(':')[1] in kegg_compounds:
                        kegg_id = id.split(':')[1]

                    # check if interactor (protein or metabolite) already exists
                    if (kegg_id is None) & (org_id != 'eco'):
                        interactor = session.query(Interactor).get(id.split(':')[1])
                        if interactor is not None:
                            interactors[num].append([interactor, None])
                    # if it doesnt exist, it's not a valid protein, so check if it is a valid compound
                    elif kegg_id is not None:
                        interactor = session.query(Metabolite).filter_by(kegg = kegg_id)
                        if interactor is None:
                            new_metabolites[num].append(kegg_id)
                        else:
                            interactors[num].append([interactor, interactor.id])
                        # if it is a valid compound, create new metabolite
                    # if parsing E. coli path, add all orthologs to interactor list
                    elif org_id == 'eco':
                        for ortholog in session.query(OrthologEcoli).filter_by(ortholog_id = id.split(':')[1],
                                                                                strain_protein = strain).all():
                            if ortholog is not None:
                                interactors[num].append([ortholog.protein, id.split(':')[1]])

            # create list of interactor pairs from two separate lists
            interactor_pairs = []
            for interactor1 in interactors[0]:
                for interactor2 in interactors[1]:
                    if (interactor1[0].type != 'm') | (interactor2[0].type != 'm'):
                        interactor_pairs.append([interactor1, interactor2])
            for interactor1 in interactors[0]:
                for id in new_metabolites[1]:
                    if interactor1[0].type == 'm': continue
                    metabolite = session.query(Metabolite).get(id)
                    if metabolite is None:
                        metabolite = Metabolite(id = id, kegg = id, pubchem = kegg_compounds[id]['pubchem'],
                                                chebi = kegg_compounds[id]['chebi'])
                        session.add(metabolite)
                    interactor_pairs.append([interactor1, [metabolite, metabolite.id]])
            for interactor1 in interactors[1]:
                for id in new_metabolites[0]:
                    if interactor1[0].type == 'm': continue
                    metabolite = session.query(Metabolite).get(id)
                    if metabolite is None:
                        metabolite = Metabolite(id = id, kegg = id, pubchem = kegg_compounds[id]['pubchem'],
                                                chebi = kegg_compounds[id]['chebi'])
                        session.add(metabolite)
                    interactor_pairs.append([interactor1, [metabolite, metabolite.id]])

            if len(interactor_pairs) == 0: continue

            # get all intermediates in reaction of type compound
            intermeds = []
            for subtype in relation.element.iter(tag='subtype'):
                if 'compound' in subtype.attrib:
                    compound_node_id = subtype.attrib['compound']
                    if compound_node_id is None: continue
                    if int(compound_node_id) not in compound_ids: continue
                    # if compound id is valid, either add existing matching metabolite or create new one and add
                    kegg_id = compound_ids[int(compound_node_id)]
                    metabolite = session.query(Metabolite).get(kegg_id)
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

            for interactor_pair in interactor_pairs:
                homogenous = (interactor_pair[0][0] == interactor_pair[1][0])
                interaction = session.query(Interaction).filter(Interaction.interactors.contains(interactor_pair[0][0]),
                                                                Interaction.interactors.contains(interactor_pair[1][0]),
                                                                Interaction.homogenous == homogenous).first()

                if interaction is None:
                    interaction = Interaction(type=interactor_pair[0][0].type + '-' + interactor_pair[1][0].type,
                                              strain=strain, homogenous=homogenous,
                                              interactors=[interactor_pair[0][0], interactor_pair[1][0]])
                    if org_id == 'eco':
                        interaction.ortholog_derived = 'Ecoli'
                    session.add(interaction), session.commit()

                interactor_a, interactor_b = None, None
                if org_id == 'eco':
                    if interaction.interactors[0] == interactor_pair[0][0]:
                        interactor_a = interactor_pair[0][1]
                        interactor_b = interactor_pair[1][1]
                    else:
                        interactor_b = interactor_pair[0][1]
                        interactor_a = interactor_pair[1][1]

                source = session.query(InteractionSource).filter_by(data_source=sourcedb).first()
                if source not in interaction.sources:
                    interaction.sources.append(source)

                reference = session.query(InteractionReference).filter_by(source_db='kegg',
                                                                          comment='in ' + path_name + ' path',
                                                                          interactor_a=interactor_a,
                                                                          interactor_b=interactor_b).first()
                if reference is None:
                    reference = InteractionReference(source_db='kegg', comment='in ' + path_name + ' path',
                                                     interactor_a=interactor_a, interactor_b=interactor_b)
                    interaction.references.append(reference)
                    reference.sources.append(source)
                else:
                    if interaction not in reference.interactions:
                        reference.interactions.append(interaction)
                    if source not in reference.sources:
                        reference.sources.append(source)

    session.commit()
    print(sourcedb, session.query(Interaction).count())