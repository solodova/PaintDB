import csv
from Schema1 import Interactor, Interaction, InteractionReference, InteractionSource, is_experimental_psimi
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

            interactor_A = session.query(Interactor).get(locus_tag1)
            if interactor_A is None: continue
            interactor_B = session.query(Interactor).get(locus_tag2)
            if interactor_B is None: continue

            homogenous = (interactors[0] == interactors[1])

            interaction = session.query(Interaction).filter(Interaction.interactors.contains(interactors[0]),
                                                            Interaction.interactors.contains(interactors[1]),
                                                             Interaction.homogenous == homogenous).first()

            if interaction is None:
                interaction = Interaction(strain='PAO1', homogenous=homogenous, interactors=interactors,
                                          type = 'p-p')
                session.add(interaction), session.commit()

            source_db, psimi_db = None, None
            if row['source_db'] != '-':
                source_db = row['source_db'].split('(')[1][:-1]
                psimi_db = row['source_db'].split('MI:')[1][:4]

            reference = session.query(InteractionReference).filter_by(psimi_detection = row['detection'].split('MI:')[1][:4],
                                                                      detection_method = row['detection'].split('(')[1][:-1],
                                                                      psimi_type = row['type'].split('MI:')[1][:4],
                                                                      interaction_type = find_type(row['type'].split(':')[2][:-1]),
                                                                      psimi_db = psimi_db,
                                                                      source_db = source_db,
                                                                      confidence = row['confidence']).first()
            if reference is None:
                reference = InteractionReference(psimi_detection=row['detection'].split('MI:')[1][:4],
                                                 detection_method=row['detection'].split('(')[1][:-1],
                                                 psimi_type=row['type'].split('MI:')[1][:4],
                                                 interaction_type=find_type(row['type'].split(':')[2][:-1]),
                                                 psimi_db=psimi_db,
                                                 source_db=source_db,
                                                 confidence=row['confidence'])
                interaction.references.append(reference)
            elif reference not in interaction.references:
                    interaction.references.append(reference)

            is_experimental = is_experimental_psimi(row['detection'].split('MI:')[1][:4])

            source = session.query(InteractionSource).filter_by(data_source = 'STRING',
                                                                is_experimental = is_experimental).first()

            if source is None:
                source = InteractionSource(data_source='STRING', is_experimental = is_experimental)
                interaction.sources.append(source)
            elif source not in interaction.sources:
                interaction.sources.append(source)


        session.commit()
