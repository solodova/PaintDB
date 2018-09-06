import csv
from Schema1 import Interactor, Interaction, InteractionReference, InteractionSource

def parse(session):
    with open('Data/PAO1/Zhang.csv') as csvfile:
        reader = csv.DictReader(csvfile)

        source = InteractionSource(data_source = 'Zhang', is_experimental = 0)
        session.add(source), session.commit()

        for row in reader:
            if float(row['Confidence']) < 0.9: continue

            interactor_A = session.query(Interactor).get(row['Protein1'])
            if interactor_A is None: continue

            interactor_B = session.query(Interactor).get(row['Protein2'])
            if interactor_B is None: continue

            homogenous = (interactor_A == interactor_B)
            interaction = session.query(Interaction).filter(Interaction.interactors.contains(interactor_A),
                                                            Interaction.interactors.contains(interactor_B),
                                                            Interaction.homogenous == homogenous).first()
            if interaction is None:
                interaction = Interaction(strain='PAO1', homogenous=homogenous, type='p-p',
                                          interactors=[interactor_A, interactor_B])
                session.add(interaction), session.commit()

            reference = session.query(InteractionReference).filter_by(detection_method = 'computational prediction',
                                                                      pmid = '22848443', interaction_type = 'predicted',
                                                                      confidence = row['Confidence'],
                                                                      comment = row['Comment']).first()
            if reference is None:
                reference = InteractionReference(detection_method='computational prediction',
                                                 author_ln='Zhang', pub_date='2012', pmid='22848443',
                                                 interaction_type='predicted', confidence=row['Confidence'],
                                                 comment=row['Comment'])
                interaction.references.append(reference)
                reference.sources.append(source)
            else:
                if reference not in interaction.references:
                    interaction.references.append(reference)
                if source not in reference.sources:
                    reference.sources.append(source)

            if source not in interaction.sources:
                interaction.sources.append(source)
    session.commit()
    print('zhang', session.query(Interaction).count())