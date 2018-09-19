import csv
from DB_schema import Interactor, Protein, Interaction, InteractionReference, InteractionSource

def parse(session):
    with open('Data/PAO1/xlinkdb.txt') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')

        reference = InteractionReference(detection_method='chemical cross-linking mass spectrometry',
                                         interaction_type='physical association', author_ln='Navari',
                                         pub_date='2015', pmid='25800553', source_db='xlinkdb')
        source = InteractionSource(data_source='XLinkDB', is_experimental = 1)
        source.references.append(reference)
        session.add(source), session.add(reference), session.commit()

        for row in reader:

            interactor_A = session.query(Interactor).get(row['proA'])
            if interactor_A is None:
                interactor_A = session.query(Protein).filter_by(uniprotkb = row['proA']).first()
            if interactor_A is None: continue

            interactor_B = session.query(Interactor).get(row['proB'])
            if interactor_B is None:
                interactor_B = session.query(Protein).filter_by(uniprotkb = row['proB']).first()
            if interactor_B is None: continue

            homogenous = (interactor_A == interactor_B)
            interaction = session.query(Interaction).filter(Interaction.interactors.contains(interactor_A),
                                                            Interaction.interactors.contains(interactor_B),
                                                            Interaction.homogenous == homogenous).first()
            if interaction is None:
                interaction = Interaction(strain='PAO1', homogenous=homogenous, type='p-p',
                                          interactors=[interactor_A, interactor_B])
                interaction.references.append(reference)
                interaction.sources.append(source)
                session.add(interaction), session.commit()
            else:
                if reference not in interaction.references:
                    interaction.references.append(reference)
                if source not in interaction.sources:
                    interaction.sources.append(source)

    session.commit()
    print('xlinkdb', session.query(Interaction).count())