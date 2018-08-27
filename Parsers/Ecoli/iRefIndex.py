import csv
from Schema1 import OrthologEcoli, Interaction, InteractionReference, InteractionSource
from Parsers.Parser import is_experimental_psimi

def parse_ecoli_irefindex(session):
    with open('Ecoli/PSICQUIC/iRefIndex.txt') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
            if (row['#ID(s) interactor A'] == '-') | (row['ID(s) interactor B'] == '-'): continue
            interactors = []

            orthologs_A = []
            id_A = row['#ID(s) interactor A'].split(':')
            if id_A[0] == 'uniprotkb':
                orthologs_A = session.query(OrthologEcoli).filter(OrthologEcoli.ortholog_uniprot == id_A[1]).all()
            elif id_A[0] == 'refseq':
                orthologs_A = session.query(OrthologEcoli).filter(OrthologEcoli.ortholog_refseq == id_A[1]).all()

            if len(orthologs_A) == 0: continue

            orthologs_B = []
            id_B = row['ID(s) interactor B'].split(':')
            if id_B[0] == 'uniprotkb':
                orthologs_B = session.query(OrthologEcoli).filter(OrthologEcoli.ortholog_uniprot == id_B[1]).all()
            elif id_B[0] == 'refseq':
                orthologs_B = session.query(OrthologEcoli).filter(OrthologEcoli.ortholog_refseq == id_B[1]).all()

            for ortholog_A in orthologs_A:
                for ortholog_B in orthologs_B:
                    if (ortholog_A is not None) and (ortholog_B is not None):
                        if ortholog_A.strain_protein == ortholog_B.strain_protein:
                            interactors.append([[ortholog_A.protein, ortholog_A.ortholog_id],
                                                [ortholog_B.protein, ortholog_B.ortholog_id]])

            for interactor_pair in interactors:
                homogenous = (interactor_pair[0][0] == interactor_pair[1][0])
                interaction = session.query(Interaction).filter(Interaction.interactors.contains(interactor_pair[0][0]),
                                                                Interaction.interactors.contains(interactor_pair[1][0]),
                                                                Interaction.homogenous == homogenous).first()
                if interaction is not None:
                    if interaction.ortholog_derived is None:
                        interaction.ortholog_derived = 'cfe'
                    elif 'fe' not in interaction.ortholog_derived:
                        interaction.ortholog_derived += ', cfe'
                    session.commit()
                else:
                    interaction = Interaction(strain=interactor_pair[0][0].strain,
                                              interactors=[interactor_pair[0][0], interactor_pair[1][0]],
                                              type=(interactor_pair[0][0].type + '-' + interactor_pair[1][0]),
                                              ortholog_derived='fe')
                    if 'MI:' in row['Interaction detection method(s)']:
                        #iterate through all methods
                        if is_experimental_psimi(row['Interaction detection method(s)'].split('MI:')[1][:4]):
                            interaction.is_experimental = 1
                    session.add(interaction), session.commit()

                interactor_a, interactor_b = '', ''
                if interaction.interactors[0] == interactor_pair[0][0]:
                    interactor_a = interactor_pair[0][1]
                    interactor_b = interactor_pair[1][1]
                else:
                    interactor_b = interactor_pair[0][1]
                    interactor_a = interactor_pair[1][1]

                author, date, psimi_type, type = None, None, None, None
                confidences, psimi_detections, detections, pmids = [None], [None], [None], [None]
                if row['Publication 1st author(s)'] != '-':
                    author = row['Publication 1st author(s)'].split(' ')[0]
                    date = row['Publication 1st author(s)'].split('(')[1][:-1]
                if row['Interaction type(s)'] != '-':
                    type = row['Interaction type(s)'].split('(')[1][:-1]
                    if 'MI' in row['Interaction type(s)']:
                        psimi_type = row['Interaction type(s)'].split('MI:')[1][:4]
                if row['Publication Identifier(s)'] != '-':
                    del pmids[0]
                    for pmid in row['Publication Identifier(s)'].split('|'):
                        pmids.append(pmid.split(':')[1])
                if row['Interaction detection method(s)'] != '-':
                    del detections[0]
                    del psimi_detections[0]
                    for detection in row['Publication Identifier(s)'].split('|'):
                        detections.append(detection.split('(')[1][:-1])
                        psimi_detections.append(detection.split('MI:')[1][:4])

                for pmid in pmids:
                    for confidence in confidences:
                        for (detection, psimi_detection) in zip(detections, psimi_detections):
                            reference = InteractionReference(interaction_id=interaction.id,
                                                             psimi_detection=psimi_detection,
                                                             detection_method=detection,
                                                             author_ln = author,
                                                             date = date,
                                                             psimi_type=psimi_type,
                                                             interaction_type=type,
                                                             psimi_db=row['Source database(s)'].split('MI')[1][:4],
                                                             source_db=row['Source database(s)'].split('(')[1][:-1],
                                                             confidence=confidence,
                                                             interactor_a=interactor_a,
                                                             interactor_b=interactor_b)
                            session.add(reference)

                source = session.query(InteractionSource).filter(
                    InteractionSource.interaction_id == interaction.id,
                    InteractionSource.data_source == 'iRefIndex').first()

                if source is None:
                    source = InteractionSource(interaction_id=interaction.id, data_source='iRefIndex')
                    session.add(source)
        session.commit()
        print(session.query(Interaction).count())