import csv
from Schema1 import OrthologEcoli, Interactor, Interaction, InteractionReference

def parse_Ecoli_RegulonDB(session):
    with open('Ecoli/RegulonDB.csv') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            interactors = []

            orthologs_A = session.query(OrthologEcoli).filter(
                OrthologEcoli.ortholog_name == (row['TF name'][0].lower() + row['TF name'][1:])).all()
            orthologs_B = session.query(OrthologEcoli).filter(
                OrthologEcoli.ortholog_name == row['Regulated gene']).all()

            for ortholog_A in orthologs_A:
                if ortholog_A != None:
                    interactor_A = session.query(Interactor).filter(Interactor.id == ortholog_A.protein_id).one()
                    for ortholog_B in orthologs_B:
                        if ortholog_B != None:
                            interactor_B = session.query(Interactor).filter(
                                Interactor.id == ortholog_B.protein_id).one()
                            if (interactor_A.strain == interactor_B.strain):
                                interactors.append([interactor_A, interactor_B])

            for interactor_pair in interactors:
                homogenous = (interactor_pair[0] == interactor_pair[1])
                interaction = session.query(Interaction).filter(
                    (Interaction.interactors.contains(interactor_pair[0])),
                    (Interaction.interactors.contains(interactor_pair[1])),
                    (Interaction.homogenous == homogenous)).first()
                if (interaction != None):
                    if (interaction.ortholog_derived == None):
                        interaction.ortholog_derived = 'confirmed from E.coli'
                    elif ('from E. coli' not in interaction.ortholog_derived):
                        interaction.ortholog_derived += ', confirmed from E. coli'
                    session.commit()
                if (interaction == None):
                    interaction = Interaction(strain=interactor_pair[0].strain, interactors=interactor_pair,
                                              type=(interactor_pair[0].type + '-' + interactor_pair[1].type),
                                              is_experimental=0, ortholog_derived='from E. coli')
                    session.add(interaction), session.commit()
                reference = InteractionReference(interaction_id=interaction.id,
                                                 interaction_type='TF/sigma-binding site (' +
                                                                  row['Regulatory effect'] + 'regulation)',
                                                 detection_method=row['Evidence'],
                                                 interaction_full_name=row['TF name'].lower() + ' regulates(' +
                                                                       row['Regulatory effect'] + ') ' +
                                                                       row['Regulated gene'],
                                                 source_db='RegulonDB')
                session.add(reference)
        session.commit()
        print(session.query(InteractionReference).count())