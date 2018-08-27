import csv
from Schema1 import Interactor, Interaction, InteractionReference, InteractionSource
from Parsers.Parser import is_experimental_psimi
def find_type(psi_code):
    return {
        '1110': 'predicted interaction',
        '2232': 'molecular association',
        '0914': 'association',
        '0915': 'physical association',
        '0407': 'direct interaction',
    }[psi_code]

def parse_string(session):
    with open('PAO1/STRING.txt') as csvfile:
        fieldnames = ['interactor_A', 'interactor_B', 'altID_A', 'altID_B', 'alias_A', 'alias_B', 'detection',
                      'publication', 'publication_ID', 'taxid_A', 'taxid_B', 'type', 'source_db', 'identifier',
                      'confidence']
        reader = csv.DictReader(csvfile, delimiter='\t', fieldnames=fieldnames)
        for row in reader:
            interactors = []
            locus_tag1 = row['interactor_A'].split('|')[0].split('.')[1]
            locus_tag2 = row['interactor_B'].split('|')[0].split('.')[1]

            if session.query(Interactor).filter(Interactor.id == locus_tag1).first() is not None:
                interactors.append(session.query(Interactor).filter(Interactor.id == locus_tag1).one())

            if session.query(Interactor).filter(Interactor.id == locus_tag2).first() is not None:
                interactors.append(session.query(Interactor).filter(Interactor.id == locus_tag2).one())

            if len(interactors) != 2: continue
            homogenous = (interactors[0] == interactors[1])

            interaction = session.query(Interaction).filter(Interaction.interactors.contains(interactors[0]),
                                                            Interaction.interactors.contains(interactors[1]),
                                                            Interaction.homogenous == homogenous).first()
            if interaction is None:
                type = interactors[0].type + '-' + interactors[1].type
                interaction = Interaction(strain='PAO1', type=type, homogenous=homogenous, interactors=interactors)
                session.add(interaction), session.commit()
                if is_experimental_psimi(row['detection'].split('MI:')[1][:4]):
                    interaction.is_experimental = 1
                else:
                    interaction.is_experimental = 0
            else:
                if is_experimental_psimi(row['detection'].split('MI:')[1][:4]):
                    interaction.is_experimental = 1

            source_db = None
            if (row['source_db'] != '-'):
                source_db = row['source_db'].split('(')[1][:-1]

            reference = InteractionReference(interaction_id=interaction.id,
                                             detection_method=row['detection'].split('(')[1][:-1],
                                             interaction_type=find_type(row['type'].split(':')[2][:-1]),
                                             source_db=source_db,
                                             confidence= row['confidence'])
            session.add(reference)

            source = session.query(InteractionSource).filter(InteractionSource.interaction_id == interaction.id,
                                                             InteractionSource.data_source == 'STRING').first()

            if source is None:
                source = InteractionSource(interaction_id=interaction.id, data_source='STRING')
                session.add(source)

        session.commit()
        print(session.query(Interaction).count())