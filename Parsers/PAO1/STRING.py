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

            interactor_A = session.query(Interactor).filter(Interactor.id == locus_tag1).first()
            interactor_B = session.query(Interactor).filter(Interactor.id == locus_tag2).first()

            if (interactor_A is None) | (interactor_B is None): continue
            homogenous = (interactors[0] == interactors[1])

            interaction = session.query(Interaction).filter(Interaction.interactors.contains(interactors[0]),
                                                            Interaction.interactors.contains(interactors[1]),
                                                             Interaction.homogenous == homogenous).first()

            is_experimental = is_experimental_psimi(row['detection'].split('MI:')[1][:4])
            if interaction is None:
                interaction = Interaction(strain='PAO1', homogenous=homogenous, interactors=interactors,
                                          type = (interactor_A.type + '-' + interactor_B.type))
                session.add(interaction), session.commit()

            source_db, psimi_db = None, None
            if row['source_db'] != '-':
                source_db = row['source_db'].split('(')[1][:-1]
                psimi_db = row['source_db'].split('MI:')[1][:4]

            reference = session.query(InteractionReference).filter(
                InteractionReference.psimi_detection == row['detection'].split('MI:')[1][:4],
                InteractionReference.detection_method == row['detection'].split('(')[1][:-1],
                InteractionReference.author_ln == None,
                InteractionReference.pub_date == None,
                InteractionReference.pmid == None,
                InteractionReference.psimi_type == row['type'].split('MI:')[1][:4],
                InteractionReference.interaction_type == find_type(row['type'].split(':')[2][:-1]),
                InteractionReference.psimi_db == psimi_db,
                InteractionReference.source_db == source_db,
                InteractionReference.confidence == row['confidence'],
                InteractionReference.comment == None,
                InteractionReference.interactor_a == None,
                InteractionReference.interactor_b == None).first()
            if reference is None:
                reference = InteractionReference(psimi_detection=row['detection'].split('MI:')[1][:4],
                                                 detection_method=row['detection'].split('(')[1][:-1],
                                                 psimi_type=row['type'].split('MI:')[1][:4],
                                                 interaction_type=find_type(row['type'].split(':')[2][:-1]),
                                                 psimi_db=psimi_db,
                                                 source_db=source_db,
                                                 confidence=row['confidence'])
                interaction.references.append(reference)
                session.add(reference), session.commit()
            else:
                if reference not in interaction.references:
                    interaction.references.append(reference)

            source = session.query(InteractionSource).filter(InteractionSource.data_source == 'STRING',
                                                             is_experimental == is_experimental).first()

            if source is None:
                source = InteractionSource(data_source='STRING', is_experimental = is_experimental)
                session.add(source), session.commit()
                interaction.sources.append(source)
            elif source not in interaction.sources:
                interaction.sources.append(source)


        session.commit()
        print(session.query(Interaction).count())