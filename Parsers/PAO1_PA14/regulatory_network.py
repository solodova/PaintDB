import csv
from Schema1 import Interactor, Protein, Interaction, InteractionReference, InteractionSource

def parse_regulatory_network(session):
    source = InteractionSource(data_source = 'Galan-Vasquez', is_experimental=2)
    session.add(source), session.commit()

    with open('Data/PAO1_PA14/regulatory_network.csv') as csvfile:
        reader = csv.DictReader(csvfile)
        string = []
        for row in reader:
            strains = row['Strain'].split(',')
            for strain in strains:
                if (strain != 'PAO1') and (strain != 'PA14'): continue
                interactors = []
                interactor_A =session.query(Protein).filter(Protein.name == row['Regulator (TF or sigma)'],
                                                            Protein.strain == strain).first()
                if interactor_A is None:
                    interactor_A=session.query(Interactor).filter(
                        Interactor.id == row['Regulator (TF or sigma)']).first()

                if interactor_A is None: continue

                interactor_B = session.query(Protein).filter(Protein.name == row['Target'],
                                                             Protein.strain == strain).first()
                if interactor_B is None:
                    interactor_B = session.query(Interactor).filter(Interactor.id == row['Target']).first()

                if interactor_B is None: continue

                string.append(len(interactors))

                # string.append(row['Regulator (TF or sigma)'] + ':' + row['Target'] + '(' + strain + ')')
                homogenous = (interactor_A == interactor_B)
                interaction = session.query(Interaction).filter(Interaction.interactors.contains(interactor_A),
                                                                Interaction.interactors.contains(interactor_B),
                                                                Interaction.homogenous == homogenous).first()

                if interaction is None:
                    interaction = Interaction(strain=strain, type='p-p', homogenous=homogenous,
                                              interactors=interactors)
                    session.add(interaction), session.commit()

                source_db, detections = None, [None]
                if row['source_db'] != '':
                    source_db=row['source_db']
                if row['detection'] != '':
                    del detections[0]
                    for type in row['detection'].split(', '):
                        detections.append(type)

                for detection in detections:
                    reference = InteractionReference(detection_method=row['evidence'],
                                                     pmid=row['pmid'],
                                                     interaction_type='TF/sigma-binding site (' + row['mode'] +
                                                                      'regulation)',
                                                     source_db=source_db,
                                                     comment=interactor_A.id + ' regulates(' +
                                                                           row['mode'] + ') ' + interactor_B.id)
                    session.add(reference)

                if source not in interaction.sources:
                    interaction.sources.append(source)

        session.commit()
        print(session.query(InteractionReference).count())