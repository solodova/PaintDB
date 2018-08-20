import csv
from Schema1 import Metabolite, Interaction, OrthologEcoli, InteractionReference, InteractionXref
from os.path import exists

ecocyc_compounds = {}

def parse_Ecoli_EcoCyc(session):
    get_ecocyc_compounds(session)
    parse_ecocyc('PAO1', session)
    parse_ecocyc('PA14', session)

def get_ecocyc_compounds(session):
    string = ''
    with open("Ecoli/ECOCYC_paths.txt") as csvfile:
        path_reader = csv.DictReader(csvfile)
        for path_row in path_reader:
            interactor_file_name = "ecocyc_files/" + path_row["Pathways"] + "_interactors.txt"
            if not exists(interactor_file_name): continue
            with open(interactor_file_name) as pathfile:
                interactor_reader = csv.DictReader(pathfile)
                for interactor_row in interactor_reader:
                    if interactor_row["PARTICIPANT_TYPE"] != "SmallMoleculeReference": continue
                    xrefs = interactor_row["UNIFICATION_XREF"]
                    name = interactor_row['PARTICIPANT']
                    if (name not in ecocyc_compounds):
                        ecocyc_compounds[name] = {}
                        id_list = [['PubChem-compound:', 'PubChem'], ['ChEBI:', 'ChEBI'], ['CAS:', 'CAS'],
                                   ['KEGG LIGAND:', 'KEGG'], ['EcoCyc:', 'EcoCyc']]
                        for id_type in id_list:
                            ecocyc_compounds[name][id_type[1]] = None
                            if (len(xrefs.split(id_type[0])) < 2): continue
                            ecocyc_compounds[name][id_type[1]] = xrefs.split(id_type[0])[1].split(';')[0]
    print(string)
    # for metabolite in ecocyc_compounds:
    #     if 'KEGG' in metabolite:
    #         if (session.query(Metabolite).filter(Metabolite.id == metabolite["KEGG"]).first() is not None):
    #             interactor = session.query(Metabolite.filter(Metabolite.id == metabolite["KEGG"])).one()
    #             interactor.CAS = metabolite['CAS']
    #             interactor.ChEBI = metabolite['ChEBI']
    #             interactor.EcoCyc = metabolite['EcoCyc']
    #             interactor.id = metabolite['EcoCyc']


def parse_ecocyc(strain, session):
    with open("Ecoli/ECOCYC_paths.txt") as csvfile:
        path_reader = csv.DictReader(csvfile)
        for path_row in path_reader:
            interaction_file_name = "ecocyc_files/" + path_row["Pathways"] + "_interactions.txt"
            if not exists(interaction_file_name): continue
            with open(interaction_file_name) as pathfile:
                interaction_reader = csv.DictReader(pathfile)
                for interaction_row in interaction_reader:
                    interactors = [[], []]
                    new_metabolites = [[], []]

                    A_id = interaction_row['PARTICIPANT_A']
                    B_id = interaction_row['PARTICIPANT_B']

                    if (A_id not in ecocyc_compounds):
                        for ortholog in (session.query(OrthologEcoli).filter(
                                (OrthologEcoli.ortholog_uniprot == A_id),
                                (OrthologEcoli.strain_protein == strain)).all()):
                            if (ortholog is not None): interactors[0].append(ortholog.protein)
                    else:
                        A_ecocyc = ecocyc_compounds[A_id]["EcoCyc"]
                        if (session.query(Metabolite).filter(Metabolite.id == A_ecocyc).first() is not None):
                            interactors[0].append(session.query(Metabolite).filter(
                                Metabolite.id == A_ecocyc).one())
                        else:
                            new_metabolites[0].append(A_ecocyc)

                    if (B_id not in ecocyc_compounds):
                        for ortholog in (session.query(OrthologEcoli).filter(
                                (OrthologEcoli.ortholog_uniprot == B_id),
                                (OrthologEcoli.strain_protein == strain)).all()):
                            if (ortholog is not None): interactors[1].append(ortholog.protein)
                    else:
                        B_ecocyc = ecocyc_compounds[B_id]["EcoCyc"]
                        if (session.query(Metabolite).filter(Metabolite.id == B_ecocyc).first() is not None):
                            interactors[1].append(session.query(Metabolite).filter(
                                Metabolite.id == B_ecocyc).one())
                        else:
                            new_metabolites[1].append(B_ecocyc)
                            print(B_id, B_ecocyc, path_row['Pathways'])

                    interactor_pairs = []
                    for interactor1 in interactors[0]:
                        for interactor2 in interactors[1]:
                            if (interactor1.type != 'metabolite') | (interactor2.type != 'metabolite'):
                                interactor_pairs.append([interactor1, interactor2])

                    for interactor1 in interactors[0]:
                        for interactor2 in new_metabolites[1]:
                            if (interactor1.type != 'metabolite'):
                                metabolite = session.query(Metabolite).filter(Metabolite.id == interactor2).first()
                                if (metabolite is None):
                                    metabolite = Metabolite(id=interactor2, EcoCyc=interactor2,
                                                            PubChem=ecocyc_compounds[B_id]['PubChem'],
                                                            KEGG=ecocyc_compounds[B_id]['KEGG'],
                                                            CAS=ecocyc_compounds[B_id]['CAS'],
                                                            ChEBI=ecocyc_compounds[B_id]['ChEBI'], name=B_id)
                                    session.add(metabolite), session.commit()
                                interactor_pairs.append([interactor1, metabolite])

                    for interactor1 in interactors[1]:
                        for interactor2 in new_metabolites[0]:
                            if (interactor1.type != 'metabolite'):
                                metabolite = session.query(Metabolite).filter(Metabolite.id == interactor2).first()
                                if (metabolite is None):
                                    metabolite = Metabolite(id=interactor2, EcoCyc=interactor2,
                                                            PubChem=ecocyc_compounds[A_id]['PubChem'],
                                                            KEGG=ecocyc_compounds[A_id]['KEGG'],
                                                            CAS=ecocyc_compounds[A_id]['CAS'],
                                                            ChEBI=ecocyc_compounds[A_id]['ChEBI'], name=A_id)
                                    session.add(metabolite), session.commit()
                                interactor_pairs.append([metabolite, interactor1])

                    if (len(interactor_pairs) == 0): continue

                    for interactor_pair in interactor_pairs:
                        homogenous = (interactor_pair[0] == interactor_pair[1])
                        interaction = session.query(Interaction).filter(
                            (Interaction.interactors.contains(interactor_pair[0])),
                            (Interaction.interactors.contains(interactor_pair[1])),
                            (Interaction.homogenous == homogenous)).first()

                        if (interaction == None):
                            interaction = Interaction(type=interactor_pair[0].type + '-' + interactor_pair[1].type,
                                                      strain=strain, homogenous=homogenous,
                                                      interactors=interactor_pair,
                                                      comment=interactor_pair[0].id +
                                                              interaction_row["INTERACTION_TYPE"] +
                                                              interactor_pair[1].id,
                                                      ortholog_derived='from E. coli')
                            session.add(interaction), session.commit()
                        else:
                            interaction = session.query(Interaction).filter(
                                (Interaction.interactors.contains(interactor_pair[0])),
                                (Interaction.interactors.contains(interactor_pair[1])),
                                (Interaction.homogenous == homogenous)).one()
                            if (interaction.ortholog_derived == None):
                                interaction.ortholog_derived = 'confirmed from E. coli'
                            else:
                                if ('from E.coli' not in interaction.ortholog_derived):
                                    interaction.ortholog_derived += ', confirmed from E.coli'

                        for pmid in interaction_row["INTERACTION_PUBMED_ID"].split(';'):
                            if (session.query(InteractionReference).filter(
                                    InteractionReference.interaction_id == interaction.id,
                                    InteractionReference.pmid == pmid).first() is None):
                                reference = InteractionReference(interaction_id=interaction.id, pmid=pmid)
                                session.add(reference), session.commit()

                        if (session.query(InteractionXref).filter(
                                (InteractionXref.interaction_id == interaction.id),
                                (InteractionXref.data_source == 'EcoCyc')).first() == None):
                            xref = InteractionXref(interaction_id=interaction.id, data_source='EcoCyc')
                            session.add(xref), session.commit()
                        session.commit()
    print(session.query(Interaction).count())