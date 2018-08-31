import csv
from Schema1 import Interactor, Interaction, InteractionSource, InteractionReference, is_experimental_psimi

def parse_geoff(session):
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
                interaction = Interaction(strain='PAO1', homogenous=homogenous,
                                          interactors = [interactor_A, interactor_B],
                                          type=(interactor_A.type + '-' + interactor_B.type))
                session.add(interaction), session.commit()

            reference = InteractionReference(detection_method=row['experimental_type'], pmid=row['pmid'])
            interaction.references.append(reference)
            session.add(reference), session.commit()

            if source not in interaction.sources:
                interaction.sources.append(source)

        session.commit()