import csv
from Schema1 import Interactor, Interaction, InteractionReference, InteractionSource

def parse_zhang(session):
    with open('Data/PAO1/Zhang.csv') as csvfile:
        reader = csv.DictReader(csvfile)

        for row in reader:
            if float(row['Confidence']) < 0.9: continue
            interactors = []
            if session.query(Interactor).filter(Interactor.id == row['Protein1']).first() is not None:
                interactors.append(session.query(Interactor).filter(Interactor.id == row['Protein1']).one())

            if session.query(Interactor).filter(Interactor.id == row['Protein2']).first() is not None:
                interactors.append(session.query(Interactor).filter(Interactor.id == row['Protein2']).one())

            if len(interactors) != 2: continue
            homogenous = (interactors[0] == interactors[1])

            interaction = session.query(Interaction).filter((Interaction.interactors.contains(interactors[0])),
                                                            (Interaction.interactors.contains(interactors[1])),
                                                            (Interaction.homogenous == homogenous)).first()
            if interaction is None:
                type = interactors[0].type + '-' + interactors[1].type
                interaction = Interaction(strain='PAO1', type=type, homogenous=homogenous, is_experimental=0,
                                          interactors=interactors)
                session.add(interaction), session.commit()

            reference = InteractionReference(interaction_id=interaction.id,
                                             detection_method='computational prediction',
                                             author_ln='Zhang',
                                             pub_date='2012',
                                             pmid = '22848443',
                                             interaction_type='predicted',
                                             confidence=row['Confidence'],
                                             comment=row['Comment'])
            session.add(reference)

            source = session.query(InteractionSource).filter(InteractionSource.interaction_id == interaction.id,
                                                             InteractionSource.data_source == 'Zhang').first()

            if source is None:
                source = InteractionSource(interaction_id=interaction.id, data_source='Zhang')
                session.add(source)

        session.commit()
    print(session.query(Interaction).count())