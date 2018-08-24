import csv
from Schema1 import Interactor, Protein, Interaction, InteractionReference, InteractionSource

def parse_xlinkdb(session):
    with open('PAO1/xlinkdb.txt') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
            interactors = []

            if session.query(Interactor).filter(Interactor.id == row['proA']).first() is not None:
                interactors.append(session.query(Interactor).filter(Interactor.id == row['proA']).one())
            elif session.query(Protein).filter(Protein.uniprotkb == row['proA']).first() is not None:
                interactors.append(session.query(Protein).filter(Protein.uniprotkb == row['proA']).one())
            if session.query(Interactor).filter(Interactor.id == row['proB']).first() is not None:
                interactors.append(session.query(Interactor).filter(Interactor.id == row['proB']).one())
            elif session.query(Protein).filter(Protein.uniprotkb == row['proB']).first() is not None:
                interactors.append(session.query(Protein).filter(Protein.uniprotkb == row['proB']).one())

            if len(interactors) != 2: continue
            homogenous = (interactors[0] == interactors[1])
            type = interactors[0].type + '-' + interactors[1].type
            alt_type = interactors[1].type + '-' + interactors[0].type
            interaction = session.query(Interaction).filter((Interaction.interactors.contains(interactors[0])),
                                                            (Interaction.interactors.contains(interactors[1])),
                                                            (Interaction.homogenous == homogenous)).first()
            if (interaction == None):
                interaction = Interaction(strain='PAO1', type=type, homogenous=homogenous, interactors=interactors,
                                          is_experimental = 1)
                session.add(interaction), session.commit()
            else:
                interaction.is_experimental = 1
                if (type not in interaction.type) and (alt_type not in interaction.type):
                    interaction.type += ', ' + type

            reference = InteractionReference(interaction_id=interaction.id,
                                             detection_method='chemical cross-linking mass spectrometry',
                                             interaction_type='physical association',
                                             pmid='25800553',
                                             source_db='XLinkDB')
            session.add(reference)

            source = session.query(InteractionSource).filter(InteractionSource.interaction_id == interaction.id,
                                                             InteractionSource.data_source == 'XLinkDB').first()
            if source is None:
                source = InteractionSource(interaction_id=interaction.id, data_source='XLinkDB')
                session.add(source)

        session.commit()
        print(session.query(Interaction).count())