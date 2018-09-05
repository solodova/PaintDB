import csv
from Schema1 import Metabolite, Interaction, OrthologEcoli, InteractionReference, InteractionSource
from os.path import exists

ecocyc_paths = []
ecocyc_compounds = {}

def parse(session):
    get_ecocyc_paths()
    get_ecocyc_compounds(session)
    update_metabolite_info(session)
    source = InteractionSource(data_source='EcoCyc', is_experimental=2)
    session.add(source), session.commit()
    parse_ecocyc('PAO1', session)
    parse_ecocyc('PA14', session)

def get_ecocyc_paths():
    with open('Data/Ecoli/ecocyc_files/paths.txt') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            ecocyc_paths.append(row['Pathways'])

def get_ecocyc_compounds(session):
    id_list = [['PubChem-compound:', 'pubchem'], ['ChEBI:', 'chebi'], ['CAS:', 'cas'], ['KEGG LIGAND:', 'kegg'],
               ['EcoCyc:', 'ecocyc']]
    for path in ecocyc_paths:
        interactor_file_name = "Data/Ecoli/ecocyc_files/interactors_sif/" + path + "_interactors.txt"
        #if there was a problem with obtaining the sif files for a pathway, they may not exist
        if not exists(interactor_file_name): continue
        with open(interactor_file_name) as file:
            reader = csv.DictReader(file)
            for row in reader:
                #ignore protein participants
                if row["PARTICIPANT_TYPE"] != "SmallMoleculeReference": continue
                xrefs = row["UNIFICATION_XREF"]
                name = row['PARTICIPANT']
                # if molecule is not in ecocyc_compounds, add all it's info
                if name not in ecocyc_compounds:
                    ecocyc_compounds[name] = {}
                    # for each metabolite id of interest, add the xref if it exists
                    for id_type in id_list:
                        ecocyc_compounds[name][id_type[1]] = None
                        if id_type[0] not in xrefs: continue
                        ecocyc_compounds[name][id_type[1]] = xrefs.split(id_type[0])[1].split(';')[0]


def update_metabolite_info(session):
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

                #if id_A isn't in ecocyc_compounds, it's probably a uniprot id; search for ecoli orthologs matching id_A
                if id_A not in ecocyc_compounds:
                    for ortholog in session.query(OrthologEcoli).filter_by(ortholog_uniprot = id_A,
                                                                           strain_protein = strain).all():
                        if ortholog is not None:
                            # add both the pseudomonas protein and the ortholog id (will be needed for reference)
                            # to interactors_A
                            interactors_A.append([ortholog.protein, ortholog.ortholog_id])
                # if id_A is in ecocyc_compounds, it means it's a metabolite id
                else:
                    A_ecocyc = ecocyc_compounds[id_A]['ecocyc']
                    #check if the metabolite already exists in database
                    metabolite = session.query(Metabolite).filter_by(ecocyc=A_ecocyc).first()
                    if metabolite is not None:
                        # if metabolite exists, add both the metabolite and it's name (needed for reference)
                        # to interactors_A
                        interactors_A.append([metabolite, metabolite.name])
                    else:
                        # if metabolite doesn't exist yet, store it's id to create it later (don't create it now
                        # since if interactor_B is invalid, there is no need for new metabolite to be created)
                        new_metabolite_A = A_ecocyc

                # same as for id_A above
                if id_B not in ecocyc_compounds:
                    for ortholog in session.query(OrthologEcoli).filter_by(ortholog_uniprot = id_B,
                                                                           strain_protein = strain).all():
                        if ortholog is not None:
                            interactors_B.append([ortholog.protein, ortholog.ortholog_id])
                else:
                    B_ecocyc = ecocyc_compounds[id_B]['ecocyc']
                    metabolite = session.query(Metabolite).filter_by(ecocyc = B_ecocyc).first()
                    if metabolite is not None:
                        # if metabolite exists, add both the metabolite and it's name (needed for reference)
                        # to interactors_B
                        interactors_B.append([metabolite, metabolite.name])
                    else:
                        new_metabolite_B = B_ecocyc

                interactor_pairs = []

                # case where no unknown metabolites were found
                if (new_metabolite_A is None) and (new_metabolite_B is None):
                    # iterate through known interactors, add them together to interactor_pairs
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
                            # and the interactor pair (for the new metabolite, make sure to add it's name (for
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
                for interactor_pair in interactor_pairs:
                    homogenous = (interactor_pair[0][0] == interactor_pair[1][0])
                    interaction = session.query(Interaction).filter(
                        Interaction.interactors.contains(interactor_pair[0][0]),
                        Interaction.interactors.contains(interactor_pair[1][0]),
                        Interaction.homogenous == homogenous).first()

                    if interaction is None:
                        interaction = Interaction(type=(interactor_pair[0][0].type + '-' + interactor_pair[1][0].type),
                                                  strain=strain, homogenous=homogenous,
                                                  interactors=[interactor_pair[0][0], interactor_pair[1][0]])
                        session.add(interaction), session.commit()

                    #in case the interaction already existed, make sure interactor_a and interactor_b variables for
                    # new interaction reference match up with the first and second interactors of the existing
                    # interaction
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
                        if reference is None:
                            reference = InteractionReference(pmid=pmid, source_db='ecocyc', comment = comment,
                                                           interactor_a=interactor_a, interactor_b=interactor_b)
                            interaction.references.append(reference)
                        elif reference not in interaction.references:
                            interaction.references.append(reference)
                    source = session.query(InteractionSource).filter_by(data_source = 'EcoCyc').first()
                    if source not in interaction.sources:
                        interaction.sources.append(source)
                session.commit()