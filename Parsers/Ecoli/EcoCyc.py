import csv
from Schema1 import Metabolite, Interaction, OrthologEcoli, InteractionReference, InteractionXref, InteractionSource
from os.path import exists
from sqlalchemy import or_

ecocyc_paths = []
ecocyc_compounds = {}

def parse_ecocyc(session):
    get_ecocyc_paths()
    get_ecocyc_compounds(session)
    parse('PAO1', session)
    #parse('PA14', session)

def get_ecocyc_paths():
    with open('Data/Ecoli/ecocyc_files/paths.txt') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            ecocyc_paths.append(row['Pathways'])

def get_ecocyc_compounds(session):
    id_list = [['PubChem-compound:', 'pubchem'], ['ChEBI:', 'chebi'], ['CAS:', 'cas'],
               ['KEGG LIGAND:', 'kegg'], ['EcoCyc:', 'ecocyc']]
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
    # for metabolite in ecocyc_compounds:
    #     if 'KEGG' in metabolite:
    #         if session.query(Metabolite).filter(Metabolite.id == metabolite["KEGG"]).first() is not None:
    #             interactor = session.query(Metabolite.filter(Metabolite.id == metabolite["KEGG"])).one()
    #             interactor.CAS = metabolite['CAS']
    #             interactor.ChEBI = metabolite['ChEBI']
    #             interactor.EcoCyc = metabolite['EcoCyc']
    #             interactor.id = metabolite['EcoCyc']


def parse(strain, session):
    source = InteractionSource(data_source = 'EcoCyc', is_experimental=2)
    session.add(source), session.commit()
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
                    for ortholog in session.query(OrthologEcoli).filter(OrthologEcoli.ortholog_uniprot == id_A,
                                                                         OrthologEcoli.strain_protein == strain).all():
                        if ortholog is not None:
                            # add both the pseudomonas protein and the ortholog id (will be needed for reference)
                            # to interactors_A
                            interactors_A.append([ortholog.protein, ortholog.ortholog_id])
                # if id_A is in ecocyc_compounds, it means it's a metabolite id
                else:
                    A_ecocyc = ecocyc_compounds[id_A]['ecocyc']
                    #check if the metabolite already exists in database
                    metabolite = session.query(Metabolite).filter(Metabolite.ecocyc == A_ecocyc).first()
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
                    for ortholog in session.query(OrthologEcoli).filter(OrthologEcoli.ortholog_uniprot == id_B,
                                                                        OrthologEcoli.strain_protein == strain).all():
                        if ortholog is not None:
                            interactors_B.append([ortholog.protein, ortholog.ortholog_id])
                else:
                    B_ecocyc = ecocyc_compounds[id_B]['ecocyc']
                    metabolite = session.query(Metabolite).filter(Metabolite.ecocyc == B_ecocyc).first()
                    if metabolite is not None:
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
                            metabolite=session.query(Metabolite).filter(Metabolite.id==new_metabolite_A).first()
                            # create a new metabolite if it doesn't exist
                            if metabolite is None:
                                metabolite = Metabolite(id=new_metabolite_A, name = id_A,
                                                        ecocyc = ecocyc_compounds[id_A]['ecocyc'],
                                                        pubchem=ecocyc_compounds[id_A]['pubchem'],
                                                        kegg=ecocyc_compounds[id_A]['kegg'],
                                                        cas= ecocyc_compounds[id_A]['cas'],
                                                        chebi=ecocyc_compounds[id_A]['chebi'])
                                session.add(metabolite), session.commit()
                            # and the interactor pair (for the new metabolite, make sure to add it's name (for
                            # reference later)
                            interactor_pairs.append([interactor_B, [metabolite, id_A]])
                # same as previous case, but if new metabolite is new_metabolite_B
                elif new_metabolite_B is not None:
                    for interactor_A in interactors_A:
                        if interactor_A[0].type != 'm':
                            metabolite=session.query(Metabolite).filter(Metabolite.id==new_metabolite_B).first()
                            if metabolite is None:
                                metabolite = Metabolite(id=new_metabolite_B, name = id_B,
                                                        ecocyc = ecocyc_compounds[id_B]['ecocyc'],
                                                        pubchem=ecocyc_compounds[id_B]['pubchem'],
                                                        kegg=ecocyc_compounds[id_B]['kegg'],
                                                        cas= ecocyc_compounds[id_B]['cas'],
                                                        chebi=ecocyc_compounds[id_B]['chebi'])
                                session.add(metabolite), session.commit()
                            interactor_pairs.append([interactor_A, [metabolite, id_B]])

                # iterate through all interactor pairs and create new interactions
                for interactor_pair in interactor_pairs:
                    homogenous = (interactor_pair[0][0] == interactor_pair[1][0])
                    interaction = session.query(Interaction).filter(
                        Interaction.interactors.contains(interactor_pair[0][0]),
                        Interaction.interactors.contains(interactor_pair[1][0]),
                        Interaction.homogenous == homogenous).first()

                    if interaction is None:
                        interaction = Interaction(type=interactor_pair[0][0].type + '-' + interactor_pair[1][0].type,
                                                  strain=strain, homogenous=homogenous, interactors=interactor_pair,
                                                  ortholog_derived='fe')
                        session.add(interaction), session.commit()
                    else:
                        # if interaction already exists but has no ortholog derivation, add 'cfe'
                        if interaction.ortholog_derived is None:
                            interaction.ortholog_derived = 'cfe'
                        # if interaction has an ortholog derivation but ecoli is not mentioned, add 'cfe'
                        elif 'fe' not in interaction.ortholog_derived:
                            interaction.ortholog_derived += ', cfe'

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

                    # iterate through all the pmids listed as reference for given interaction
                    for pmid in interaction_row["INTERACTION_PUBMED_ID"].split(';'):
                        comment = interactor_pair[0][1] + interaction_row["INTERACTION_TYPE"] + interactor_pair[1][1]
                        # check if interaction reference already exists in db
                        reference = session.query(InteractionReference).filter(
                                InteractionReference.psimi_detection == None,
                                InteractionReference.detection_method == None,
                                InteractionReference.author_ln == None,
                                InteractionReference.pub_date == None,
                                InteractionReference.pmid == pmid,
                                InteractionReference.psimi_type == None,
                                InteractionReference.interaction_type == None,
                                InteractionReference.psimi_db == None,
                                InteractionReference.source_db == 'ecocyc',
                                InteractionReference.confidence == None,
                                InteractionReference.comment == comment,
                                InteractionReference.interactor_a == interactor_a,
                                InteractionReference.interactor_b == interactor_b).first()
                        if reference is None:
                            new_ref = InteractionReference(pmid=pmid, source_db='ecocyc', comment = comment,
                                                           interactor_a=interactor_a, interactor_b=interactor_b)
                            session.add(new_ref), session.commit()
                            interaction.references.append(new_ref)
                        elif reference not in interaction.references:
                            interaction.references.append(reference)

                    if source not in interaction.sources:
                        interaction.sources.append(source)
                session.commit()
    print(session.query(Interaction).count())
    #print(session.query(Interaction).filter(Interaction.type == 'p-p').count())