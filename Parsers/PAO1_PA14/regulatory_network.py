import csv
from Schema import Interactor, Protein, Interaction, InteractionReference, InteractionSource

def parse(session):
    # create and add sources for the interactions (do this before since they all use the same source)
    # Note: is_experimental is set to 2 because we cannot confirm that detection method was experimental or not
    source_PAO1 = InteractionSource(data_source = 'Galan-Vasquez(PAO1)', is_experimental=2)
    source_PA14 = InteractionSource(data_source='Galan-Vasquez(PA14)', is_experimental=2)
    session.add(source_PAO1), session.add(source_PA14), session.commit()

    with open('Data/PAO1_PA14/regulatory_network.csv') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            # a row describing an interaction may have >1 strain
            strains = row['Strain'].split(',')
            for strain in strains:
                # only care about PAO1 and PA14 strain interactions
                if (strain != 'PAO1') and (strain != 'PA14'): continue

                # search for interactor A by name
                interactor_A =session.query(Protein).filter_by(name = row['Regulator'], strain = strain).first()
                # if no interactor was found by name, id listed may be a gene locus, so search by this id
                if interactor_A is None:
                    interactor_A=session.query(Interactor).get(row['Regulator'])
                # if no interactor A was found for this interaction, skip to next
                if interactor_A is None: continue

                # same as A above
                interactor_B = session.query(Protein).filter_by(name = row['Target'], strain = strain).first()
                if interactor_B is None:
                    interactor_B = session.query(Interactor).get(row['Target'])
                if interactor_B is None: continue

                homogenous = (interactor_A == interactor_B)
                interaction = session.query(Interaction).filter(Interaction.interactors.contains(interactor_A),
                                                                Interaction.interactors.contains(interactor_B),
                                                                Interaction.homogenous == homogenous).first()
                # if interaction between these 2 interactors does not yet exist, create and add it
                if interaction is None:
                    interaction = Interaction(strain=strain, type='p-p', homogenous=homogenous,
                                              interactors=[interactor_A, interactor_B])
                    session.add(interaction), session.commit()

                # specify the source to be used for the interaction and reference based on strain of interaction
                source = None
                if strain == 'PAO1':
                    source = source_PAO1
                else:
                    source = source_PA14

                # add the source to the interaction source list if it isn't there already
                if source not in interaction.sources:
                    interaction.sources.append(source)

                # get source db and detections if they are present in the file
                source_db, detections = None, [None]
                if row['source_db'] != '':
                    source_db = row['source_db']
                if row['evidence'] != '':
                    del detections[0]
                    for type in row['evidence'].split(', '):
                        detections.append(type)
                # create a new reference for each detection found, add the reference to the interaction's
                # reference list, and add the source to the reference's sources
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