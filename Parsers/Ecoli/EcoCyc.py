import csv
from Schema import Metabolite, Interaction, OrthologEcoli, InteractionReference, InteractionSource
from os.path import exists

# list of ecocyc paths
ecocyc_paths = []
# dict of compounds in ecocyc interactor files, based on interactors found in all pathways
# {name: {pubchem: None, chebi: None, cas: None, kegg: None, ecocyc: None}}
ecocyc_compounds = {}

# main function to go through process of parsing EcoCyc
def parse(session):
    # parse ecocyc paths from file from EcoCyc and pyut into ecocyc_paths
    get_ecocyc_paths()
    # parse compounds from ecocyc interactor files, put into ecocyc_compounds
    get_ecocyc_compounds(session)
    update_metabolite_info_ecocyc(session)
    # create and add new source for EcoCyc (no references for any interactions so is_experimental = 2)
    source = InteractionSource(data_source='EcoCyc', is_experimental=2)
    session.add(source), session.commit()
    # parse PAO1 and PA14 separately
    parse_ecocyc('PAO1', session)
    parse_ecocyc('PA14', session)

# function to get ecocyc paths from file
def get_ecocyc_paths():
    with open('Data/Ecoli/ecocyc_files/paths.txt') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            ecocyc_paths.append(row['Pathways'])

# function to parse compounds from ecocyc interactor files, put into ecocyc_compounds
def get_ecocyc_compounds(session):
    id_list = [['PubChem-compound:', 'pubchem'], ['ChEBI:', 'chebi'], ['CAS:', 'cas'], ['KEGG LIGAND:', 'kegg'],
               ['EcoCyc:', 'ecocyc']]
    for path in ecocyc_paths:
        interactor_file_name = "Data/Ecoli/ecocyc_files/interactors_sif/" + path + "_interactors.txt"
        # if there was a problem with obtaining the interactor sif files for a pathway, they may not exist
        if not exists(interactor_file_name): continue
        # if file does exist, open it and parse
        with open(interactor_file_name) as file:
            reader = csv.DictReader(file)
            for row in reader:
                # ignore protein participants, only want to collect compounds for compound dict
                if row["PARTICIPANT_TYPE"] != "SmallMoleculeReference": continue
                xrefs = row["UNIFICATION_XREF"]
                name = row['PARTICIPANT']
                # if molecule is not in ecocyc_compounds, add all it's info
                if name not in ecocyc_compounds:
                    ecocyc_compounds[name] = {}
                    # for each metabolite id of interest (in id_type), add the xref if it exists
                    for id_type in id_list:
                        # set id to None in case it doesn't exist
                        ecocyc_compounds[name][id_type[1]] = None
                        if id_type[0] not in xrefs: continue
                        # if the id is present, add it to ecocyc dict[name]
                        ecocyc_compounds[name][id_type[1]] = xrefs.split(id_type[0])[1].split(';')[0]


# function to go through each compound in ecocyc_compounds and try to find it in PaIntDB using each id present in dict,
# if it exists then update all missing attributes
# note: by calling this function before parsing EcoCyc, existing metabolites will be assigned ecocyc ids
# (if existing metabolite matches any ids in ecocyc_compounds); thus metabolites can be searched by their
# ecocyc id alone in parse_ecocyc
def update_metabolite_info_ecocyc(session):
    for name in ecocyc_compounds:
        metabolite = None
        kegg, pubchem, ecocyc, chebi, cas =  ecocyc_compounds[name]['kegg'],  ecocyc_compounds[name]['pubchem'], \
                                             ecocyc_compounds[name]['ecocyc'],  ecocyc_compounds[name]['chebi'],  \
                                             ecocyc_compounds[name]['cas']
        if pubchem is not None:
            metabolite = session.query(Metabolite).filter_by(pubchem = pubchem).first()
        if (metabolite is None) & (chebi is not None):
            metabolite = session.query(Metabolite).filter_by(chebi=chebi).first()
        if (metabolite is None) & (kegg is not None):
            metabolite = session.query(Metabolite).filter_by(kegg=kegg).first()
        if (metabolite is None) & (cas is not None):
            metabolite = session.query(Metabolite).filter_by(cas=cas).first()
        if (metabolite is None) & (ecocyc is not None):
            metabolite = session.query(Metabolite).filter_by(ecocyc=ecocyc).first()
        if metabolite is None:
            metabolite = session.query(Metabolite).filter_by(name=name).first()

        if metabolite is not None:
            if metabolite.kegg is None:
                metabolite.kegg = kegg
            if metabolite.pubchem is None:
                metabolite.pubchem = pubchem
            if metabolite.chebi is None:
                metabolite.chebi = chebi
            if metabolite.cas is None:
                metabolite.cas = cas
            if metabolite.ecocyc is None:
                metabolite.ecocyc = ecocyc
            if metabolite.name is None:
                metabolite.name = name


# function to parse interactions from EcoCyc with ortholog mapping to given strain (PAO1 or PA14)
def parse_ecocyc(strain, session):
    for path in ecocyc_paths:
        interaction_file_name = "Data/Ecoli/ecocyc_files/interactions_sif/" + path + "_interactions.txt"
        #if there was a problem with obtaining the sif files for a pathway, they may not exist
        if not exists(interaction_file_name): continue
        with open(interaction_file_name) as file:
            reader = csv.DictReader(file)
            for interaction_row in reader:
                interactors_A, interactors_B = [], []
                new_metabolite_A, new_metabolite_B = None, None

                id_A = interaction_row['PARTICIPANT_A']
                id_B = interaction_row['PARTICIPANT_B']

                #if id_A isn't in ecocyc_compounds, it's a uniprot id; search for ecoli orthologs matching id_A
                if id_A not in ecocyc_compounds:
                    for ortholog in session.query(OrthologEcoli).filter_by(ortholog_uniprot = id_A,
                                                                           strain_protein = strain).all():
                        if ortholog is not None:
                            # add both the pseudomonas protein and the ortholog id (will be needed later to
                            # create interaction reference) to interactors_A
                            interactors_A.append([ortholog.protein, ortholog.ortholog_id])
                # if id_A is in ecocyc_compounds, it means it's a metabolite id
                else:
                    A_ecocyc = ecocyc_compounds[id_A]['ecocyc']
                    #check if the metabolite already exists in database (only need to search ecocyc id since
                    # update_metabolites_ecocyc was called)
                    metabolite = session.query(Metabolite).filter_by(ecocyc=A_ecocyc).first()
                    if metabolite is not None:
                        # if metabolite exists, add both the metabolite and it's name (will be needed later
                        # to create interaction reference) to interactors_A
                        interactors_A.append([metabolite, metabolite.name])
                    else:
                        # if metabolite doesn't exist yet, store it's id to create it later (don't create it now
                        # since if interactor_B is invalid, there is no need for new metabolite to be created)
                        new_metabolite_A = A_ecocyc

                # same as for id_A above, now with second interactor
                if id_B not in ecocyc_compounds:
                    for ortholog in session.query(OrthologEcoli).filter_by(ortholog_uniprot = id_B,
                                                                           strain_protein = strain).all():
                        if ortholog is not None:
                            interactors_B.append([ortholog.protein, ortholog.ortholog_id])
                else:
                    B_ecocyc = ecocyc_compounds[id_B]['ecocyc']
                    metabolite = session.query(Metabolite).filter_by(ecocyc = B_ecocyc).first()
                    if metabolite is not None:
                        interactors_B.append([metabolite, metabolite.name])
                    else:
                        new_metabolite_B = B_ecocyc

                # store new interactor pairs from which to create interactions here
                interactor_pairs = []

                # case where no unknown metabolites were found
                if (new_metabolite_A is None) and (new_metabolite_B is None):
                    # iterate through known protein interactors, add them together to interactor_pairs
                    for interactor_A in interactors_A:
                        for interactor_B in interactors_B:
                            # only add the interactor pair if at least one of them is not a metabolite
                            if (interactor_A[0].type != 'm') | (interactor_B[0].type != 'm'):
                                interactor_pairs.append([interactor_A, interactor_B])
                # case where there is one new metabolite (new_metabolite_A)
                elif new_metabolite_A is not None:
                    for interactor_B in interactors_B:
                        # don't add a new interactor pair if both are metabolites
                        if interactor_B[0].type != 'm':
                            # check if new metabolite exists in database (eg. if more than one ortholog was found for
                            # interactors_B, you don't want to create the same new metabolite twice)
                            metabolite=session.query(Metabolite).filter_by(ecocyc=new_metabolite_A).first()
                            # create a new metabolite if it doesn't exist
                            if metabolite is None:
                                metabolite = Metabolite(id=new_metabolite_A, name = id_A,
                                                        ecocyc = new_metabolite_A,
                                                        pubchem=ecocyc_compounds[id_A]['pubchem'],
                                                        kegg=ecocyc_compounds[id_A]['kegg'],
                                                        cas= ecocyc_compounds[id_A]['cas'],
                                                        chebi=ecocyc_compounds[id_A]['chebi'])
                                session.add(metabolite)
                            # add the interactor pair (for the new metabolite, make sure to add it's name (for
                            # reference later)
                            interactor_pairs.append([interactor_B, [metabolite, id_A]])
                # same as previous case, but if new metabolite is new_metabolite_B
                elif new_metabolite_B is not None:
                    for interactor_A in interactors_A:
                        if interactor_A[0].type != 'm':
                            metabolite=session.query(Metabolite).filter_by(ecocyc=new_metabolite_B).first()
                            if metabolite is None:
                                metabolite = Metabolite(id=new_metabolite_B, name = id_B,
                                                        ecocyc = new_metabolite_B,
                                                        pubchem=ecocyc_compounds[id_B]['pubchem'],
                                                        kegg=ecocyc_compounds[id_B]['kegg'],
                                                        cas= ecocyc_compounds[id_B]['cas'],
                                                        chebi=ecocyc_compounds[id_B]['chebi'])
                                session.add(metabolite)
                            interactor_pairs.append([interactor_A, [metabolite, id_B]])

                # iterate through all interactor pairs and create new interactions
                # note interactor_pairs will be empty if:
                #   1) both interactors were new metabolites
                #   2) one or both ecoli interactors did not have orthologs in Pseudomonas
                for interactor_pair in interactor_pairs:
                    homogenous = (interactor_pair[0][0] == interactor_pair[1][0])
                    interaction = session.query(Interaction).filter(
                        Interaction.interactors.contains(interactor_pair[0][0]),
                        Interaction.interactors.contains(interactor_pair[1][0]),
                        Interaction.homogenous == homogenous).first()

                    source = session.query(InteractionSource).filter_by(data_source='EcoCyc').first()

                    # if interaction doesn't exist, add it, and EcoCyc as a source
                    if interaction is None:
                        # if this interaction is created for first time, mark it as ortholog derived from Ecoli
                        interaction = Interaction(type=(interactor_pair[0][0].type + '-' + interactor_pair[1][0].type),
                                                  strain=strain, homogenous=homogenous,
                                                  interactors=[interactor_pair[0][0], interactor_pair[1][0]],
                                                  ortholog_derived = 'Ecoli')
                        interaction.sources.append(source)
                        session.add(interaction), session.commit()
                    # add EcoCyc as source for interaction if it isn't already
                    elif source not in interaction.sources:
                        interaction.sources.append(source)

                    # in case the interaction already existed, make sure interactor_a and interactor_b variables for
                    # new interaction reference match up with the first and second interactors of the existing
                    # interaction (so it's easy to see which Pseudomonas interactor matches up with which Ecoli
                    # ortholog)
                    interactor_a, interactor_b = None, None
                    if interaction.interactors[0] == interactor_pair[0][0]:
                        interactor_a = interactor_pair[0][1]
                        interactor_b = interactor_pair[1][1]
                    else:
                        interactor_b = interactor_pair[0][1]
                        interactor_a = interactor_pair[1][1]

                    comment = interactor_pair[0][1] + interaction_row["INTERACTION_TYPE"] + interactor_pair[1][1]

                    # iterate through all the pmids listed as reference for given interaction
                    for pmid in interaction_row["INTERACTION_PUBMED_ID"].split(';'):
                        # check if interaction reference already exists in db
                        reference = session.query(InteractionReference).filter_by(pmid = pmid, source_db = 'ecocyc',
                                                                                  comment = comment,
                                                                                  interactor_a = interactor_a,
                                                                                  interactor_b = interactor_b).first()
                        # if reference doesn't exist, create it, add the interaction to its references, and the
                        # EcoCyc source to its sources
                        if reference is None:
                            reference = InteractionReference(pmid=pmid, source_db='ecocyc', comment = comment,
                                                           interactor_a=interactor_a, interactor_b=interactor_b)
                            reference.interactions.append(interaction)
                            reference.sources.append(source)
                        # if reference does exist, add interaction to its interactions and source to its sources
                        # (if it doesn't have them already)
                        else:
                            if interaction not in reference.interactions:
                                reference.interactions.append(interaction)
                            if source not in reference.sources:
                                reference.sources.append(source)
    session.commit()
    print('ecocyc', session.query(Interaction).count())
