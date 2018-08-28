import csv
from Schema1 import Interactor, Interaction, InteractionSource, InteractionReference, is_experimental_psimi

def parse_geoff(session):
    with open('Data/PAO1/GeoffWinsor.csv') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            interactors = []

            if session.query(Interactor).filter(Interactor.id == row['locus_tag']).first() is not None:
                interactors.append(session.query(Interactor).filter(Interactor.id == row['locus_tag']).one())
            next(reader)
            if session.query(Interactor).filter(Interactor.id == row['locus_tag']).first() is not None:
                interactors.append(session.query(Interactor).filter(Interactor.id == row['locus_tag']).one())

            if len(interactors) != 2: continue
            homogenous = (interactors[0] == interactors[1])

            interaction = session.query(Interaction).filter(Interaction.interactors.contains(interactors[0]),
                                                            Interaction.interactors.contains(interactors[1]),
                                                            Interaction.homogenous == homogenous).first()
            if interaction is None:
                type = interactors[0].type + '-' + interactors[1].type
                interaction = Interaction(strain='PAO1', type=type, homogenous=homogenous, interactors=interactors)
                session.add(interaction), session.commit()


            reference = InteractionReference(detection_method=row['experimental type'],
                                             pmid=row['pmid'])
            if session.query(InteractionReference).filter(InteractionReference == reference).first() is None:
                interaction.references.append(reference)
                session.add(reference), session.commit()
            else:
                if reference not in interaction.references:
                    interaction.references.append(reference)

            source = session.query(InteractionSource).filter(InteractionSource.interaction_id == interaction.id,
                                                             InteractionSource.data_source == 'Geoff').first()

            if source is None:
                source = InteractionSource(interaction_id=interaction.id, data_source='Geoff', is_experimental=1)
                session.add(source)

        session.commit()
        print(session.query(Interaction).count())