import csv
from Schema1 import Interactor, Interaction, InteractionReference, InteractionSource

def parse_zhang(session):
    with open('Data/PAO1/Zhang.csv') as csvfile:
        reader = csv.DictReader(csvfile)

        source = InteractionSource(data_source = 'Zhang', is_experimental = 0)
        session.add(source), session.commit()

        for row in reader:
            if float(row['Confidence']) < 0.9: continue

            interactor_A = session.query(Interactor).filter(Interactor.id == row['Protein1']).first()
            if interactor_A is None: continue

            interactor_B = session.query(Interactor).filter(Interactor.id == row['Protein2']).first()
            if interactor_B is None: continue

            homogenous = (interactor_A == interactor_B)

            interaction = session.query(Interaction).filter(Interaction.interactors.contains(interactor_A),
                                                            Interaction.interactors.contains(interactor_B),
                                                            Interaction.homogenous == homogenous).first()
            if interaction is None:
                interaction = Interaction(strain='PAO1', homogenous=homogenous, is_experimental=0,
                                          interactors=[interactor_A, interactor_B],
                                          type=(interactor_A.type + '-' + interactor_B.type))
                session.add(interaction), session.commit()

            reference = session.query(InteractionReference).filter(
                InteractionReference.psimi_detection == None,
                InteractionReference.detection_method == 'computational prediction',
                InteractionReference.author_ln == 'Zhang',
                InteractionReference.pub_date == '2012',
                InteractionReference.pmid == '22848443',
                InteractionReference.psimi_type == None,
                InteractionReference.interaction_type == 'predicted',
                InteractionReference.psimi_db == None,
                InteractionReference.source_db == None,
                InteractionReference.confidence == row['Confidence'],
                InteractionReference.comment == row['Comment'],
                InteractionReference.interactor_a == None,
                InteractionReference.interactor_b == None).first()
            if reference is None:
                reference = InteractionReference(detection_method='computational prediction',
                                                 author_ln='Zhang',
                                                 pub_date='2012',
                                                 pmid='22848443',
                                                 interaction_type='predicted',
                                                 confidence=row['Confidence'],
                                                 comment=row['Comment'])
                session.add(reference), session.commit()
                interaction.references.append(reference)
            elif reference not in interaction.references:
                interaction.references.append(reference)

            if source not in interaction.sources:
                interaction.sources.append(source)

        session.commit()
    print(session.query(Interaction).count())