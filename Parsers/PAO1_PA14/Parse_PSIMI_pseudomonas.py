import csv
from Schema1 import Interactor, Protein, Interaction, InteractionReference, InteractionSource, InteractionXref
from Schema1 import is_experimental_psimi

import itertools

cols = ['interactor_A', 'interactor_B', 'altID_A', 'altID_B', 'alias_A', 'alias_B', 'detection', 'publication',
        'publication_ID', 'taxid_A', 'taxid_B', 'type', 'source_db', 'identifier', 'confidence']

def parse(session):
    parse_psimi('Data/PAO1/PSICQUIC/ADIPInteractomes.txt', 'PAO1', 'ADIPInteractomes(PAO1)', session)
    parse_psimi('Data/PAO1/PSICQUIC/IMEx.txt', 'PAO1', 'IMEx(PAO1)', session)
    parse_psimi('Data/PAO1/PSICQUIC/IntAct.txt', 'PAO1', 'IntAct(PAO1)', session)
    parse_psimi('Data/PAO1/PSICQUIC/iRefIndex.txt', 'PAO1', 'iRefIndex(PAO1)', session)
    parse_psimi('Data/PAO1/PSICQUIC/mentha.txt', 'PAO1', 'mentha(PAO1)', session)
    parse_psimi('Data/PAO1/PSICQUIC/MINT.txt', 'PAO1', 'MINT(PAO1)', session)
    parse_psimi('Data/PAO1/PSICQUIC/MPIDB.txt', 'PAO1', 'MPIDB(PAO1)', session)
    parse_psimi('Data/PAO1/IntAct.txt', 'PAO1', 'IntAct(PAO1)', session)

    parse_psimi('Data/PA14/PSICQUIC/IMEx.txt', 'PA14', 'IMEx(PA14)', session)
    parse_psimi('Data/PA14/PSICQUIC/iRefIndex.txt', 'PA14', 'iRefIndex(PA14)', session)
    parse_psimi('Data/PA14/PSICQUIC/mentha.txt', 'PA14', 'mentha(PA14)', session)
    parse_psimi('Data/PA14/IntAct.txt', 'PA14', 'IntAct(PA14)', session)

def parse_psimi(file, strain, source, session):
    with open(file) as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t', fieldnames=cols)
        next(reader)
        for row in reader:
            uniprot_A, refseq_A, interactor_A, uniprot_B, refseq_B, interactor_B = None, None, None, None, None, None

            if 'uniprotkb' in row['interactor_A']:
                uniprot_A = row['interactor_A'].split('uniprotkb:')[1].split('|')[0]
            if 'refseq' in row['interactor_A']:
                refseq_A = row['interactor_A'].split('refseq:')[1].split('|')[0]

            if uniprot_A is not None:
                interactor_A = session.query(Interactor).get(uniprot_A)
                if interactor_A is None:
                    interactor_A = session.query(Protein).filter_by(uniprotkb = uniprot_A).first()
            if (interactor_A is None) and (refseq_A is not None):
                interactor_A = session.query(Protein).filter_by(ncbi_acc = refseq_A).first()
            if interactor_A is None: continue

            if 'uniprotkb' in row['interactor_B']:
                uniprot_B = row['interactor_B'].split('uniprotkb:')[1].split('|')[0]
            if 'refseq' in row['interactor_B']:
                refseq_B = row['interactor_B'].split('refseq:')[1].split('|')[0]

            if uniprot_B is not None:
                interactor_B = session.query(Interactor).get(uniprot_B)
                if interactor_B is None:
                    interactor_B = session.query(Protein).filter_by(uniprotkb = uniprot_B).first()
            if (interactor_B is None) and (refseq_B is not None):
                interactor_B = session.query(Protein).filter_by(ncbi_acc = refseq_B).first()

            if interactor_B is None: continue
            homogenous = (interactor_A == interactor_B)

            interaction = session.query(Interaction).filter(Interaction.interactors.contains(interactor_A),
                                                            Interaction.interactors.contains(interactor_B),
                                                            Interaction.homogenous == homogenous).first()
            if interaction is None:
                interaction = Interaction(strain=strain, type='p-p', homogenous=homogenous,
                                          interactors=[interactor_A, interactor_B])
                session.add(interaction), session.commit()

            ref_fields = {'detections': [], 'types': [], 'dbs': [], 'confidences': [], 'authors': [], 'dates': [],
                          'pmids': [], 'psimi_detections': [], 'psimi_types': [], 'psimi_dbs': []}

            if 'MI' in row['detection']:
                for psimi_detection in row['detection'].split('MI:')[1:]:
                    if psimi_detection == '': continue
                    ref_fields['psimi_detections'].append(psimi_detection[:4])
            if 'MI' in row['type']:
                for psimi_type in row['type'].split('MI:')[1:]:
                    if psimi_type == '': continue
                    ref_fields['psimi_types'].append(psimi_type[:4])
            if 'MI' in row['source_db']:
                for psimi_db in row['source_db'].split('MI:')[1:]:
                    if psimi_db == '': continue
                    ref_fields['psimi_dbs'].append(psimi_db[:4])

            for detection in row['detection'].split('|'):
                if (detection == '-') | (detection == ''): continue
                ref_fields['detections'].append(detection.split('(')[1][:-1])
            for pub in row['publication'].split('|'):
                if (pub == '-') | (pub == ''): continue
                seps = [' ', '(']
                author = ''
                if ('-' in pub) and (' ' not in pub):
                    seps = ['-', '-']
                if '_' in pub:
                    author = pub.split(seps[0])[0][0].upper() + pub.split(seps[0])[0].split('_')[0][1:]+ '-' + \
                             pub.split(seps[0])[0].split('_')[1]
                else:
                    author = pub.split(seps[0])[0][0].upper() + pub.split(seps[0])[0][1:]
                ref_fields['authors'].append(author)
                if (seps[1] == '-'):
                    ref_fields['dates'].append(pub.split(seps[1])[1])
                elif ('(' in pub):
                    ref_fields['dates'].append(pub.split(seps[1])[1][:-1])
            for id in row['publication_ID'].split('|'):
                if ('pubmed' not in id) | (id == '-') | (id == '') | ('DIP' in id): continue
                ref_fields['pmids'].append(id.split('pubmed:')[1])
            for type in row['type'].split('|'):
                if (type == '-') | (type == ''): continue
                ref_fields['types'].append(type.split('(')[1][:-1])
            for db in row['source_db'].split('|'):
                if (db == '-') | (db == ''): continue
                ref_fields['dbs'].append(db.split('(')[1][:-1])
            for confidence in row['confidence'].split('|'):
                if (confidence == '-') | (confidence == ''): continue
                if (confidence.split(':')[0] == 'core') | (confidence.split(':')[0] == 'ist'): continue
                ref_fields['confidences'].append(confidence)

            for field in ref_fields:
                if len(ref_fields[field]) == 0:
                    ref_fields[field].append(None)

            if (ref_fields['psimi_detections'][0] is None) and (ref_fields['detections'][0] is not None):
                for i in range(1, len(ref_fields['detections'])):
                    ref_fields['psimi_detections'].append(None)
            if (ref_fields['psimi_types'][0] is None) and (ref_fields['types'][0] is not None):
                for i in range(1, len(ref_fields['types'])):
                    ref_fields['psimi_types'].append(None)
            if (ref_fields['psimi_dbs'][0] is None) and (ref_fields['dbs'][0] is not None):
                for i in range(1, len(ref_fields['dbs'])):
                    ref_fields['psimi_dbs'].append(None)

            detections_full = list(zip(ref_fields['psimi_detections'], ref_fields['detections']))
            types_full = list(zip(ref_fields['psimi_types'], ref_fields['types']))
            dbs_full = list(zip(ref_fields['psimi_dbs'], ref_fields['dbs']))

            if (ref_fields['dates'][0] is None) and (ref_fields['authors'][0] is not None):
                for i in range(1, len(ref_fields['authors'])):
                    ref_fields['dates'].append(None)

            if (ref_fields['authors'][0] is None) and (ref_fields['dates'][0] is None) and \
                    (ref_fields['pmids'][0] is not None):
                for i in range(1, len(ref_fields['pmids'])):
                    ref_fields['dates'].append(None)
                    ref_fields['authors'].append(None)
            if (ref_fields['authors'][0] is not None) and (ref_fields['dates'][0] is not None) and \
                    (ref_fields['pmids'][0] is None):
                for i in range(1, len(ref_fields['authors'])):
                    ref_fields['pmids'].append(None)
            if (ref_fields['authors'][0] is not None) and (ref_fields['dates'][0] is None):
                if ref_fields['pmids'][0] is None:
                    for i in range(1, len(ref_fields['authors'])):
                        ref_fields['dates'].append(None)
                        ref_fields['pmids'].append(None)
                else:
                    if len(ref_fields['authors']) == len(ref_fields['pmids']):
                        ref_fields['dates'].append(None)

            pub_full = [(None, None, None)]
            if (len(ref_fields['authors']) == len(ref_fields['dates'])) and \
                    (len(ref_fields['authors']) == len(ref_fields['pmids'])):
                pub_full = list(zip(ref_fields['authors'], ref_fields['dates'], ref_fields['pmids']))
            if (len(list(detections_full)) > 1) & (len(list(pub_full)) == 1) & (len(list(types_full)) == 1):
                for i in range(1, len(detections_full)):
                    pub_full.append(pub_full[0])
                    types_full.append(types_full[0])

            ref_full = list(zip(detections_full, pub_full, types_full))

            ref_parameter_list = []

            for comb in itertools.product(ref_full, dbs_full, ref_fields['confidences']):
                ref_parameters = [comb[0][0][0], comb[0][0][1], comb[0][1][0], comb[0][1][1], comb[0][1][2],
                                  comb[0][2][0], comb[0][2][1], comb[1][0], comb[1][1], comb[2]]
                if all(parameter is None for parameter in ref_parameters): continue
                ref_parameter_list.append(ref_parameters)

            for ref in ref_parameter_list:
                nref = session.query(InteractionReference).filter_by(psimi_detection = ref[0],
                                                                     detection_method = ref[1], author_ln = ref[2],
                                                                     pub_date = ref[3], pmid = ref[4],
                                                                     psimi_type = ref[5], interaction_type = ref[6],
                                                                     psimi_db = ref[7], source_db = ref[8],
                                                                     confidence = ref[9]).first()
                if nref is None:
                    nref = InteractionReference(psimi_detection=ref[0], detection_method=ref[1], author_ln=ref[2],
                                                     pub_date=ref[3], pmid=ref[4], psimi_type=ref[5],
                                                     interaction_type=ref[6], psimi_db=ref[7], source_db=ref[8],
                                                     confidence=ref[9])
                    interaction.references.append(nref)
                elif nref not in interaction.references:
                    interaction.references.append(nref)

            is_experimental, not_experimental, experimental = 2, 2, 2
            for psimi_detection in ref_fields['psimi_detections']:
                if psimi_detection is not None:
                    if is_experimental_psimi(psimi_detection):
                        is_experimental = 1
                    else:
                        not_experimental = 1
            if is_experimental == 1:
                experimental = 1
            elif not_experimental == 1:
                experimental = 0

            new_source = session.query(InteractionSource).filter_by(data_source = source,
                                                                    is_experimental = experimental).first()

            if new_source is None:
                new_source = InteractionSource(data_source=source, is_experimental=experimental)
                interaction.sources.append(new_source)
            elif new_source not in interaction.sources:
                interaction.sources.append(new_source)

            for xref in row['identifier'].split('|'):
                xref_field = xref.split(':')
                xref = session.query(InteractionXref).filter_by(accession = xref_field[1],
                                                                interaction_id = interaction.id).first()

                if xref is None:
                    xref = InteractionXref(interaction_id=interaction.id, accession=xref_field[1],
                                           data_source=xref_field[0])
                    session.add(xref)

        session.commit()