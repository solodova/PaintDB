import csv
from DB_schema import Interactor, Interaction, InteractionSource, InteractionReference

def parse(session):
    with open('Data/PAO1/GeoffWinsor.csv') as csvfile:
        reader = csv.DictReader(csvfile)

        source = InteractionSource(data_source='Geoff', is_experimental=1)
        session.add(source), session.commit()

        for row in reader:
            interactor_A = session.query(Interactor).get(row['locus_tag'])
            if interactor_A is None: continue
            row = next(reader)
            interactor_B = session.query(Interactor).get(row['locus_tag'])
            if interactor_B is None: continue

            homogenous = (interactor_A == interactor_B)
            interaction = session.query(Interaction).filter(Interaction.interactors.contains(interactor_A),
                                                            Interaction.interactors.contains(interactor_B),
                                                            Interaction.homogenous == homogenous).first()
            if interaction is None:
                interaction = Interaction(strain='PAO1', homogenous=homogenous , type='p-p',
                                          interactors = [interactor_A, interactor_B])
                interaction.sources.append(source)
                session.add(interaction), session.commit()
            elif source not in interaction.sources:
                interaction.sources.append(source)

            reference = session.query(InteractionReference).filter_by(detection_method=row['experimental_type'],
                                                                      pmid=row['pmid']).first()

            if reference is None:
                reference = InteractionReference(detection_method=row['experimental_type'], pmid=row['pmid'])
                interaction.references.append(reference)
                reference.sources.append(source)
            else:
                if interaction not in reference.interactions:
                    reference.interactions.append(interaction)
                if source not in reference.sources:
                    reference.sources.append(source)

    session.commit()
    print('geoff', session.query(Interaction).count())