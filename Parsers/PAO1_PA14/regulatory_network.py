import csv
from Schema import Interactor, Protein, Interaction, InteractionReference, InteractionSource

def parse(session):
    source_PAO1 = InteractionSource(data_source = 'Galan-Vasquez(PAO1)', is_experimental=2)
    source_PA14 = InteractionSource(data_source='Galan-Vasquez(PA14)', is_experimental=2)
    session.add(source_PAO1), session.add(source_PA14), session.commit()

    with open('Data/PAO1_PA14/regulatory_network.csv') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            strains = row['Strain'].split(',')
            for strain in strains:
                if (strain != 'PAO1') and (strain != 'PA14'): continue

                interactor_A =session.query(Protein).filter_by(name = row['Regulator'], strain = strain).first()
                if interactor_A is None:
                    interactor_A=session.query(Interactor).get(row['Regulator'])
                if interactor_A is None: continue

                interactor_B = session.query(Protein).filter_by(name = row['Target'], strain = strain).first()
                if interactor_B is None:
                    interactor_B = session.query(Interactor).get(row['Target'])
                if interactor_B is None: continue

                # string.append(row['Regulator (TF or sigma)'] + ':' + row['Target'] + '(' + strain + ')')
                homogenous = (interactor_A == interactor_B)
                interaction = session.query(Interaction).filter(Interaction.interactors.contains(interactor_A),
                                                                Interaction.interactors.contains(interactor_B),
                                                                Interaction.homogenous == homogenous).first()

                if interaction is None:
                    interaction = Interaction(strain=strain, type='p-p', homogenous=homogenous,
                                              interactors=[interactor_A, interactor_B])
                    session.add(interaction), session.commit()

                source = None
                if strain == 'PAO1':
                    source = source_PAO1
                else:
                    source = source_PA14

                if source not in interaction.sources:
                    interaction.sources.append(source)

                source_db, detections = None, [None]
                if row['source_db'] != '':
                    source_db = row['source_db']
                if row['evidence'] != '':
                    del detections[0]
                    for type in row['evidence'].split(', '):
                        detections.append(type)
                for detection in detections:
                    reference = InteractionReference(detection_method=detection, pmid=row['pmid'],
                                                     interaction_type='TF/sigma-binding site (' + row['mode'] +
                                                                      'regulation)',
                                                     source_db=source_db,
                                                     comment=interactor_A.id + ' regulates(' +
                                                                           row['mode'] + ') ' + interactor_B.id)
                    interaction.references.append(reference)
                    reference.sources.append(source)

    session.commit()
    print('regnet', session.query(Interaction).count())