import csv
from Schema1 import Interactor, Metabolite, Interaction, InteractionReference, OrthologEcoli, InteractionXref, \
    InteractionSource
from Schema1 import is_experimental_psimi
import itertools


cols = ['interactor_A', 'interactor_B', 'altID_A', 'altID_B', 'alias_A', 'alias_B', 'detection', 'publication',
        'publication_ID', 'taxid_A', 'taxid_B', 'type', 'source_db', 'identifier', 'confidence']

def parse(session):
    parse_psimi(session, 'Data/Ecoli/PSICQUIC/BindingDB.txt', 'BindingDB(Ecoli)')
    parse_psimi(session, 'Data/Ecoli/PSICQUIC/EBI-GOA-nonIntAct.txt', 'EBI-GOA-nonIntAct(Ecoli)')
    parse_psimi(session, 'Data/Ecoli/PSICQUIC/IMEx.txt', 'IMEx(Ecoli)')
    parse_psimi(session, 'Data/Ecoli/PSICQUIC/IntAct.txt', 'IntAct(Ecoli)')
    parse_psimi(session, 'Data/Ecoli/PSICQUIC/iRefIndex.txt', 'iRefIndex(Ecoli)')
    parse_psimi(session, 'Data/Ecoli/PSICQUIC/mentha.txt', 'mentha(Ecoli)')
    parse_psimi(session, 'Data/Ecoli/PSICQUIC/MINT.txt', 'MINT(Ecoli)')
    parse_psimi(session, 'Data/Ecoli/PSICQUIC/MPIDB.txt', 'MPIDB(Ecoli)')
    parse_psimi(session, 'Data/Ecoli/PSICQUIC/UniProt.txt', 'UniProt(Ecoli)')
    parse_psimi(session, 'Data/Ecoli/DIP.txt', 'DIP(Ecoli)')
    parse_psimi(session, 'Data/Ecoli/IntAct.txt', 'IntAct(Ecoli)')

def parse_psimi(session, file, source):
    with open(file) as csvfile:
        reader = csv.DictReader(csvfile, fieldnames=cols, delimiter = '\t')

        # iterate through each interaction
        for row in reader:
            #if (row['interactor_A'] == '-') | (row['interactor_B'] == '-'): continue
            uniprot_A, refseq_A, orthologs_A , uniprot_B, refseq_B, orthologs_B= None, None, None, None, None, None
            pubchem, chebi = None, None
            metabolite_info, metabolite, orthologs = None, None, None

            if 'uniprotkb' in row['interactor_A']:
                uniprot_A = row['interactor_A'].split('uniprotkb:')[1].split('|')[0]
            if 'refseq' in row['interactor_A']:
                refseq_A = row['interactor_A'].split('refseq:')[1].split('|')[0]

            if uniprot_A is not None:
                orthologs_A = session.query(OrthologEcoli).filter_by(ortholog_uniprot = uniprot_A).all()
            if (orthologs_A is None) and (refseq_A is not None):
                orthologs_A = session.query(OrthologEcoli).filter_by(ortholog_refseq = refseq_A).all()
            if (orthologs_A is None) & ((uniprot_A is not None) | (refseq_A is not None)): continue

            if 'uniprotkb' in row['interactor_B']:
                uniprot_B = row['interactor_B'].split('uniprotkb:')[1].split('|')[0]
            if 'refseq' in row['interactor_B']:
                refseq_B = row['interactor_B'].split('refseq:')[1].split('|')[0]

            if uniprot_B is not None:
                orthologs_B = session.query(OrthologEcoli).filter_by(ortholog_uniprot = uniprot_B).all()
            if (orthologs_B is None) and (refseq_B is not None):
                orthologs_B = session.query(OrthologEcoli).filter_by(ortholog_refseq = refseq_B).all()
            if (orthologs_B is None) & ((uniprot_B is not None) | (refseq_B is not None)): continue

            if (orthologs_A is None) and (orthologs_B is None): continue

            if orthologs_A is None:
                if 'chebi' in row['interactor_A']:
                    chebi = row['interactor_A'].split('CHEBI:')[1].split('|')[0][:-1]
                if 'pubchem' in row['altID_A']:
                    pubchem = row['altID_A'].split('pubchem:')[1].split('|')[0]
                if (chebi is None) & ('chebi' in row['altID_A']):
                    chebi = row['altID_A'].split('CHEBI:')[1].split('|')[0][:-1]
                if (chebi is None) & (pubchem is None): continue
                orthologs = orthologs_B
            elif orthologs_B is None:
                if 'chebi' in row['interactor_B']:
                    chebi = row['interactor_B'].split('CHEBI:')[1].split('|')[0][:-1]
                if 'pubchem' in row['altID_B']:
                    pubchem = row['altID_B'].split('pubchem:')[1].split('|')[0]
                if (chebi is None) & ('chebi' in row['altID_B']):
                    chebi = row['altID_B'].split('CHEBI:')[1].split('|')[0][:-1]
                if (chebi is None) & (pubchem is None): continue
                orthologs = orthologs_A

            if (chebi is not None) | (pubchem is not None):
                id = None
                if chebi is not None:
                    id = chebi
                    metabolite = session.query(Metabolite).filter_by(chebi = chebi).first()
                if (metabolite is None) & (pubchem is not None):
                    id = pubchem
                    metabolite = session.query(Metabolite).filter_by(pubchem = pubchem).first()
                if metabolite is None:
                    metabolite = Metabolite(id = id, chebi=chebi, pubchem=pubchem)
                    session.add(metabolite)

            interactors = []
            if metabolite is None:
                for ortholog_A in orthologs_A:
                    for ortholog_B in orthologs_B:
                        if (ortholog_A is not None) and (ortholog_B is not None):
                            if ortholog_A.strain_protein == ortholog_B.strain_protein:
                                interactors.append([[ortholog_A.protein, ortholog_A.ortholog_id],
                                                    [ortholog_B.protein, ortholog_B.ortholog_id]])
            else:
                for ortholog in orthologs:
                    interactors.append([[metabolite, metabolite.id], [ortholog.protein, ortholog.ortholog_id]])

            for interactor_pair in interactors:
                homogenous = (interactor_pair[0][0] == interactor_pair[1][0])
                interaction = session.query(Interaction).filter(Interaction.interactors.contains(interactor_pair[0][0]),
                                                                Interaction.interactors.contains(interactor_pair[1][0]),
                                                                Interaction.homogenous == homogenous).first()
                if interaction is None:
                    strain = None
                    if interactor_pair[0][0].type == 'p':
                        strain = interactor_pair[0][0].strain
                    else:
                        strain = interactor_pair[1][0].strain
                    interaction = Interaction(strain=strain,
                                              interactors=[interactor_pair[0][0], interactor_pair[1][0]],
                                              type=(interactor_pair[0][0].type + '-' + interactor_pair[1][0].type),
                                              ortholog_derived='Ecoli')
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
                        author = pub.split(seps[0])[0][0].upper() + pub.split(seps[0])[0].split('_')[0][1:] + '-' + \
                                 pub.split(seps[0])[0].split('_')[1]
                    else:
                        author = pub.split(seps[0])[0][0].upper() + pub.split(seps[0])[0][1:]
                    ref_fields['authors'].append(author)
                    if seps[1] == '-':
                        ref_fields['dates'].append(pub.split(seps[1])[1])
                    elif '(' in pub:
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

                interactor_a, interactor_b = None, None
                if interaction.interactors[0] == interactor_pair[0][0]:
                    interactor_a = interactor_pair[0][1]
                    interactor_b = interactor_pair[1][1]
                else:
                    interactor_b = interactor_pair[0][1]
                    interactor_a = interactor_pair[1][1]

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

                nsource = session.query(InteractionSource).filter_by(
                    data_source=source, is_experimental=experimental).first()
                if nsource is None:
                    nsource = InteractionSource(data_source=source, is_experimental=experimental)
                    interaction.sources.append(nsource)
                elif nsource not in interaction.sources:
                    interaction.sources.append(nsource)

                for ref in ref_parameter_list:
                    nref = session.query(InteractionReference).filter_by(
                        psimi_detection=ref[0],detection_method=ref[1], author_ln = ref[2], pub_date=ref[3],
                        pmid=ref[4], psimi_type=ref[5], interaction_type=ref[6], psimi_db=ref[7], source_db=ref[8],
                        confidence=ref[9], interactor_a=interactor_a, interactor_b=interactor_b).first()
                    if nref is None:
                        nref = InteractionReference(
                            psimi_detection=ref[0], detection_method=ref[1], author_ln=ref[2], pub_date=ref[3],
                            pmid=ref[4], psimi_type=ref[5], interaction_type=ref[6], psimi_db=ref[7], source_db=ref[8],
                            confidence=ref[9], interactor_a=interactor_a, interactor_b=interactor_b)
                        interaction.references.append(nref)
                        nref.sources.append(nsource)
                    else:
                        if interaction not in nref.interactions:
                            nref.interactions.append(interaction)
                        if nsource not in nref.sources:
                            nref.sources.append(nsource)

    session.commit()
    print(source, session.query(Interaction).count())