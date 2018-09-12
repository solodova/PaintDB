import csv
from Schema import Interactor, Protein, Interaction, InteractionReference, InteractionSource, InteractionXref
from Schema import is_experimental_psimi

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

                # dict to hold all fields for references
            ref_fields = {'detections': [], 'types': [], 'dbs': [], 'confidences': [], 'authors': [], 'dates': [],
                          'pmids': []}

            # collect detection method(s), publication info, interaction type(s), source db(s)
            for detection in row['detection'].split('|'):
                if (detection == '-') | (detection == ''): continue
                ref_fields['detections'].append(detection)
            for pub in row['publication'].split('|'):
                if (pub == '-') | (pub == ''): continue
                # default separators for reference of form Zhang et al. (2014)
                # reference
                seps = [' ', '(']
                author = ''
                # if a '-' is in the reference and there are no spaces, the reference is of form Zhang_li-2014-1
                if ('-' in pub) and (' ' not in pub):
                    seps = ['-', '-']
                if '_' in pub:
                    # '_' separates combined last names in references using '-' sep
                    author = pub.split(seps[0])[0][0].upper() + pub.split(seps[0])[0].split('_')[0][1:] + '-' + \
                             pub.split(seps[0])[0].split('_')[1]
                else:
                    author = pub.split(seps[0])[0][0].upper() + pub.split(seps[0])[0][1:]
                if author[-1] == ',':
                    author = author[:-1]
                ref_fields['authors'].append(author)
                # add date based on sep
                if seps[1] == '-':
                    ref_fields['dates'].append(pub.split(seps[1])[1])
                elif '(' in pub:
                    ref_fields['dates'].append(pub.split(seps[1])[1][:-1])
            for id in row['publication_ID'].split('|'):
                # only care about pubmed ids
                if ('pubmed' not in id) | (id == '-') | (id == '') | ('DIP' in id): continue
                ref_fields['pmids'].append(id.split('pubmed:')[1])
            for type in row['type'].split('|'):
                if (type == '-') | (type == ''): continue
                ref_fields['types'].append(type)
            for db in row['source_db'].split('|'):
                if (db == '-') | (db == ''): continue
                ref_fields['dbs'].append(db)
            for confidence in row['confidence'].split('|'):
                if (confidence == '-') | (confidence == ''): continue
                # don't care about these confidence types
                if (confidence.split(':')[0] == 'core') | (confidence.split(':')[0] == 'ist'): continue
                ref_fields['confidences'].append(confidence)

            # if no value(s) were found for field, set first item of field list in ref_fields to None
            for field in ref_fields:
                if len(ref_fields[field]) == 0:
                    ref_fields[field].append(None)

            # the lengths of the author, pmid, and date lists in ref_fields should all match
            # if one or more of these lists have no values, their lengths should be extended with a None value
            # to match the length of the not-empty list

            # if there are pmids but no authors or dates, extend authors and dates
            if (ref_fields['authors'][0] is None) and (ref_fields['dates'][0] is None) and \
                    (ref_fields['pmids'][0] is not None):
                for i in range(1, len(ref_fields['pmids'])):
                    ref_fields['dates'].append(None)
                    ref_fields['authors'].append(None)
            # if there are authors and dates but no pmids, extend pmids
            elif (ref_fields['authors'][0] is not None) and (ref_fields['dates'][0] is not None) and \
                    (ref_fields['pmids'][0] is None):
                for i in range(1, len(ref_fields['authors'])):
                    ref_fields['pmids'].append(None)
            # if there are authors but no dates or pmids, extend dates and pmids
            elif (ref_fields['authors'][0] is not None) and (ref_fields['dates'][0] is None) and \
                    (ref_fields['pmids'][0] is None):
                for i in range(1, len(ref_fields['authors'])):
                    ref_fields['dates'].append(None)
                    ref_fields['pmids'].append(None)
            # if there are authors and pmids but no dates, extend dates
            elif (ref_fields['authors'][0] is not None) and (ref_fields['dates'][0] is None) and \
                    (ref_fields['pmids'][0] is not None):
                for i in range(1, len(ref_fields['authors'])):
                    ref_fields['dates'].append(None)

            # list of tuples of author, date, pmid
            pub_full = [(None, None, None)]
            # only add full info for publication if the lengths of author, date, and pmid lists match
            # Note: above list extentions should make the lists the same length; if the lists all had
            # non-None values, then this will check that all authors and dates match to a pmid;
            # if they don't then pub_full will remain empty since no 1:1:1 mapping was possible
            if (len(ref_fields['authors']) == len(ref_fields['dates'])) and \
                    (len(ref_fields['authors']) == len(ref_fields['pmids'])):
                pub_full = list(zip(ref_fields['authors'], ref_fields['dates'], ref_fields['pmids']))
            # the cases for ratio of detection to publication to type are: n:1:1, 1:1:1, n:n:n.
            # When ratio is n:1:1, the publication and interaction types lists must be extended
            # so all lengths match
            if (len(ref_fields['detections']) > 1) & (len(pub_full) == 1) & (len(ref_fields['types']) == 1):
                for i in range(1, len(ref_fields['detections'])):
                    pub_full.append(pub_full[0])
                    ref_fields['types'].append(ref_fields['types'][0])

            # create list of full references, where each reference has a detection, publication, interaction type
            ref_full = list(zip(ref_fields['detections'], pub_full, ref_fields['types']))

            # here store full list of reference parameters, including confidence and dbs
            ref_parameter_list = []

            # generate all possible cross-product combinations of the full reference (detection, publication,
            # interaction type), databases, and confidences
            # Note: this is done because the confidence scores and databases listed for an interaction
            # apply to all the full references present for the interaction
            for comb in itertools.product(ref_full, ref_fields['dbs'], ref_fields['confidences']):
                # comb looks like: (detection, (author, date, pmid), type), db, confidence
                ref_parameters = [comb[0][0], comb[0][1][0], comb[0][1][1], comb[0][1][2], comb[0][2], comb[1], comb[2]]
                # dont add the ref_parametrs to the ref_parameter list if all parameters are None
                if all(parameter is None for parameter in ref_parameters): continue
                ref_parameter_list.append(ref_parameters)

            # go through all psimi detection codes (if there are any) and determine if they are experimental
            is_experimental, not_experimental, experimental = 2, 2, 2
            if 'MI' in row['detections']:
                for psimi_detection in row['detections'].split('MI:')[1:]:
                    if psimi_detection == '': continue
                    if is_experimental_psimi(psimi_detection[:4]):
                        experimental = 1
                    else:
                        not_experimental = 1

            # if at least one detection code was experimental, experimental will be 1; set the source's
            # is_experimental attribute to 1
            if experimental == 1:
                is_experimental = 1
            # if experimental was not flagged, but not_experimental was, then set source's is_experimental
            # attribute to 0
            # Note: if both experimental and non-experimental = 1, then there were both types of detection codes
            # for this interaction, but experimental will trump non-experimental in terms of the source
            elif not_experimental == 1:
                is_experimental = 0

            # check to see if source exists
            nsource = session.query(InteractionSource).filter_by(
                data_source=source, is_experimental=is_experimental).first()
            # if source doesn't exist, create and add it to the interaction's sources
            if nsource is None:
                nsource = InteractionSource(data_source=source, is_experimental=is_experimental)
                interaction.sources.append(nsource)
            # if the source does exist, add it to the interaction's sources if it isn't already
            elif nsource not in interaction.sources:
                interaction.sources.append(nsource)

            # go through each reference in the ref_parameter list, search for it, and if it doesnt exist create it
            for ref in ref_parameter_list:
                nref = session.query(InteractionReference).filter_by(
                    detection_method=ref[0], author_ln=ref[1], pub_date=ref[2], pmid=ref[3],
                    interaction_type=ref[4], source_db=ref[5], confidence=ref[6],
                    interactor_a=None, interactor_b=None).first()
                # if nref doesn't exist, create and add it to the interaction's reference list,
                # and add the source to the reference's sources
                if nref is None:
                    nref = InteractionReference(
                        detection_method=ref[0], author_ln=ref[1], pub_date=ref[2], pmid=ref[3],
                        interaction_type=ref[4], source_db=ref[5], confidence=ref[6])
                    interaction.references.append(nref)
                    nref.sources.append(nsource)
                # if nref does exist, add the interaction and source to it's attributes if they aren't added
                else:
                    if interaction not in nref.interactions:
                        nref.interactions.append(interaction)
                    if nsource not in nref.sources:
                        nref.sources.append(nsource)

            for xref in row['identifier'].split('|'):
                xref_field = xref.split(':')
                xref = session.query(InteractionXref).filter_by(accession = xref_field[1],
                                                                interaction_id = interaction.id).first()

                if xref is None:
                    xref = InteractionXref(interaction_id=interaction.id, accession=xref_field[1],
                                           data_source=xref_field[0])
                    session.add(xref)

        session.commit()
    print(source, session.query(Interaction).count())