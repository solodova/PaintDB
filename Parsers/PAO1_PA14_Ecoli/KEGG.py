from Bio.KEGG.REST import kegg_list, kegg_get, kegg_conv
from Bio.KEGG.KGML.KGML_parser import read
from Schema1 import Interactor, Metabolite, OrthologEcoli, Interaction, InteractionReference, InteractionSource

kegg_compounds = {}

def parse_pseudomonas_KEGG(session):
    parse_KEGG('pae', 'PAO1', session)
    parse_KEGG('pau', 'PA14', session)

def parse_Ecoli_KEGG(session):
    parse_KEGG('eco', 'PAO1', session)
    parse_KEGG('eco', 'PA14', session)

def find_type_KEGG(attrib):
    return {
        'ECrel': 'enzyme-enzyme relation',
        'PPrel': 'protein-protein interaction',
        'GErel': 'gene expression interaction',
        'PCrel': 'protein-compound interaction',
    }[attrib]

# make sure to run this before calling parse_?_KEGG()
def get_KEGG_compounds():
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

def parse_KEGG(org_id, strain, session):
    # get pathways for organism specified by org_id
    pathways = kegg_list(database='pathway', org=org_id).read().split('path:')
    path_names, path_ids = [], []

    for path in pathways:
        if path != '':
            path_names.append(path.split('\t')[1].split(' -')[0])
            path_ids.append(path[:8])

    source = InteractionSource(data_source = 'KEGG', is_experimental = 2)
    session.add(source), session.commit()

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
                        interactor = session.query(Interactor).filter(Interactor.id == id.split(':')[1]).first()
                        if interactor is not None:
                            interactors[num].append([interactor])
                    # if it doesnt exist, it's not a valid protein, so check if it is a valid compound
                    elif kegg_id is not None:
                        interactor = session.query(Metabolite).filter(Metabolite.kegg == kegg_id).first()
                        if interactor is None:
                            new_metabolites[num].append(kegg_id)
                        else:
                            interactors[num].append([interactor, interactor.name])
                        # if it is a valid compound, create new metabolite
                    # if parsing E. coli path, add all orthologs to interactor list
                    elif org_id == 'eco':
                        for ortholog in (session.query(OrthologEcoli).filter(
                                OrthologEcoli.ortholog_id == id.split(':')[1],
                                OrthologEcoli.strain_protein == strain).all()):
                            interactors[num].append([ortholog.protein, ortholog.ortholog_id])

            # create list of interactor pairs from two separate lists (interactors_sif[0], interactors_sif[1])
            interactor_pairs = []
            for interactor1 in interactors[0]:
                for interactor2 in interactors[1]:
                    if (interactor1[0].type != 'm') | (interactor2[0].type != 'm'):
                        interactor_pairs.append([interactor1, interactor2])
            for interactor1 in interactors[0]:
                for id in new_metabolites[1]:
                    if interactor1[0].type == 'm': continue
                    metabolite = session.query(Metabolite).filter(Metabolite.kegg == id).first()
                    if metabolite is None:
                        metabolite = Metabolite(id = id, kegg = id, pubchem = kegg_compounds[id]['pubchem'],
                                                chebi = kegg_compounds[id]['chebi'])
                        session.add(metabolite), session.commit()
                    interactor_pairs.append([interactor1, [metabolite, metabolite.name]])
            for interactor1 in interactors[1]:
                for id in new_metabolites[0]:
                    if interactor1[0].type == 'm': continue
                    metabolite = session.query(Metabolite).filter(Metabolite.kegg == id).first()
                    if metabolite is None:
                        metabolite = Metabolite(id = id, kegg = id, pubchem = kegg_compounds[id]['pubchem'],
                                                chebi = kegg_compounds[id]['chebi'])
                        session.add(metabolite), session.commit()
                    interactor_pairs.append([interactor1, [metabolite, metabolite.name]])

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
                    metabolite = session.query(Metabolite).filter(Metabolite.kegg == kegg_id).first()
                    if metabolite is None:
                        metabolite = Metabolite(id=kegg_id, name=kegg_compounds[kegg_id]['name'],
                                                pubchem=kegg_compounds[kegg_id]['pubchem'],
                                                chebi=kegg_compounds[kegg_id]['chebi'])
                        session.add(metabolite), session.commit()
                    intermeds.append([metabolite, metabolite.name])

            # add protein - intermediate interactor pairs
            for interactor_list in interactors:
                for interactor in interactor_list[0]:
                    if interactor.type != 'm':
                        for intermed in intermeds:
                            interactor_pairs.append([interactor, [intermed, intermed.name]])

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
                        interaction.ortholog_derived = 'fe'

                    session.add(interaction), session.commit()
                if org_id == 'eco':
                    if interaction.ortholog_derived is None:
                        interaction.ortholog_derived = 'cfe'
                    elif 'fe' not in interaction.ortholog_derived:
                        interaction.ortholog_derived += ', cfe'

                    interactor_a, interactor_b = None, None
                    if interaction.interactors[0] == interactor_pair[0][0]:
                        interactor_a = interactor_pair[0][1]
                        interactor_b = interactor_pair[1][1]
                    else:
                        interactor_b = interactor_pair[0][1]
                        interactor_a = interactor_pair[1][1]

                    reference = session.query(InteractionReference).filter(
                        InteractionReference.psimi_detection == None,
                        InteractionReference.detection_method == None,
                        InteractionReference.author_ln == None,
                        InteractionReference.pub_date == None,
                        InteractionReference.pmid == None,
                        InteractionReference.psimi_type == None,
                        InteractionReference.interaction_type == None,
                        InteractionReference.psimi_db == None,
                        InteractionReference.source_db == 'kegg',
                        InteractionReference.confidence == None,
                        InteractionReference.comment == None,
                        InteractionReference.interactor_a == interactor_a,
                        InteractionReference.interactor_b == interactor_b).first()
                    if reference is None:
                        reference = InteractionReference(source_db = 'kegg',
                                                         interactor_a=interactor_a, interactor_b=interactor_b)
                        session.add(reference), session.commit()
                        interaction.references.append(reference)
                    elif reference not in interaction.references:
                        interaction.references.append(reference)

                if source not in interaction.sources:
                    interaction.sources.append(source)

    print(session.query(Interaction).count())