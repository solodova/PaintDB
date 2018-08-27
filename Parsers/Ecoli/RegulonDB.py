import csv
from Schema1 import OrthologEcoli, Interactor, Interaction, InteractionReference, InteractionSource

def parse_ecoli_regulondb(session):
    with open('Ecoli/RegulonDB.csv') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            interactors = []

            orthologs_A = session.query(OrthologEcoli).filter(
                OrthologEcoli.ortholog_name == (row['TF name'][0].lower() + row['TF name'][1:])).all()
            orthologs_B = session.query(OrthologEcoli).filter(
                OrthologEcoli.ortholog_name == row['Regulated gene']).all()

            for ortholog_A in orthologs_A:
                for ortholog_B in orthologs_B:
                    if (ortholog_A is not None) and (ortholog_B is not None):
                        if (ortholog_A.strain_protein == ortholog_B.strain_protein):
                            interactors.append([[ortholog_A.protein, ortholog_A.ortholog_id],
                                                [ortholog_B.protein, ortholog_B.ortholog_id]])

            for interactor_pair in interactors:
                homogenous = (interactor_pair[0][0] == interactor_pair[1][0])
                interaction = session.query(Interaction).filter(Interaction.interactors.contains(interactor_pair[0][0]),
                                                                Interaction.interactors.contains(interactor_pair[1][0]),
                                                                Interaction.homogenous == homogenous).first()
                if interaction is not None:
                    if interaction.ortholog_derived is None:
                        interaction.ortholog_derived = 'cfe'
                    elif 'fe' not in interaction.ortholog_derived:
                        interaction.ortholog_derived += ', cfe'
                    session.commit()
                else:
                    interaction = Interaction(strain=interactor_pair[0][0].strain,
                                              interactors=[interactor_pair[0][0], interactor_pair[1][0]],
                                              type=(interactor_pair[0][0].type + '-' + interactor_pair[1][0]),
                                              ortholog_derived='fe')

                for evidence in row['Evidence'][1:-1].split(', '):
                    reference = InteractionReference(interaction_id=interaction.id,
                                                     detection_method=evidence,
                                                     interaction_type='TF/sigma-binding site (' +
                                                                      row['Regulatory effect'] + 'regulation)',
                                                     confidence=row['Evidence type'],
                                                     comment=interactor_pair[0][1] + ' regulates(' +
                                                                           row['Regulatory effect'] + ') ' +
                                                                           interactor_pair[1][1],
                                                     source_db='RegulonDB')
                    session.add(reference)

                source = session.query(InteractionSource).filter(InteractionSource.interaction_id == interaction.id,
                                                                 InteractionSource.data_source == 'RegulonDB').first()

                if source is None:
                    source = InteractionSource(interaction_id=interaction.id, data_source='RegulonDB')
                    session.add(source)

        session.commit()
        print(session.query(InteractionReference).count())