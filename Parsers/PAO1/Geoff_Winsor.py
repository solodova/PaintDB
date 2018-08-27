import csv
from Schema1 import Interactor, Interaction, InteractionSource, InteractionReference

def parse_geoff(session):
    with open('Data/PAO1/GeoffWinsor.csv') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            interactors = []

            if session.query(Interactor).filter(Interactor.id == row['locus_tag']).first() is not None:
                interactors.append(session.query(Interactor).filter(Interactor.id == row['locus_tag']).one())
            row = next(reader)
            if session.query(Interactor).filter(Interactor.id == row['locus_tag']).first() is not None:
                interactors.append(session.query(Interactor).filter(Interactor.id == row['locus_tag']).one())

            if len(interactors) != 2: continue
            homogenous = (interactors[0] == interactors[1])

            interaction = session.query(Interaction).filter(Interaction.interactors.contains(interactors[0]),
                                                            Interaction.interactors.contains(interactors[1]),
                                                            Interaction.homogenous == homogenous).first()
            if interaction is None:
                type = interactors[0].type + '-' + interactors[1].type
                interaction = Interaction(strain='PAO1', type=type, homogenous=homogenous, interactors=interactors,
                                          is_experimental=1)
                session.add(interaction), session.commit()
            else:
                interaction.is_experimental = 1

            reference = InteractionReference(interaction_id=interaction.id,
                                             detection_method=row['experimental type'],
                                             pmid=row['pmid'],
                                             comment=row['full_name'])
            session.add(reference)

            source = session.query(InteractionSource).filter(InteractionSource.interaction_id == interaction.id,
                                                             InteractionSource.data_source == 'Geoff').first()

            if source is None:
                source = InteractionSource(interaction_id=interaction.id, data_source='Geoff')
                session.add(source)

        session.commit()
        print(session.query(Interaction).count())