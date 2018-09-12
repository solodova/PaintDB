import csv
from Schema import Interactor, Metabolite, Interaction, InteractionReference, OrthologEcoli, InteractionXref, \
    InteractionSource
from Schema import is_experimental_psimi
import itertools

# universal column names to be used for all file parsing
cols = ['interactor_A', 'interactor_B', 'altID_A', 'altID_B', 'alias_A', 'alias_B', 'detection', 'publication',
        'publication_ID', 'taxid_A', 'taxid_B', 'type', 'source_db', 'identifier', 'confidence']

# function to parse all Ecoli psimi formatted data files
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

            uniprot_A, refseq_A, orthologs_A , uniprot_B, refseq_B, orthologs_B= None, None, None, None, None, None
            # if one of the interactors is metabolte, save it's ids in pubchem and chebi
            pubchem, chebi = None, None
            # if one of the interactors is a metabolite, metabolite will be that metabolite and orthologs
            # will be set to the interaction's protein ortholog(s)
            metabolite_info, metabolite, orthologs = None, None, None

            # check if interactor A has uniprot or refseq id
            if 'uniprotkb' in row['interactor_A']:
                uniprot_A = row['interactor_A'].split('uniprotkb:')[1].split('|')[0]
            if 'refseq' in row['interactor_A']:
                refseq_A = row['interactor_A'].split('refseq:')[1].split('|')[0]

            # if uniprot id was found, look for orthologs matching that id
            if uniprot_A is not None:
                orthologs_A = session.query(OrthologEcoli).filter_by(ortholog_uniprot = uniprot_A).all()
            # if no orthologs were found but a refseq id was found, try to find ortholog based on refseq
            if (orthologs_A is None) and (refseq_A is not None):
                orthologs_A = session.query(OrthologEcoli).filter_by(ortholog_refseq = refseq_A).all()
            # if no orthologs were found for interactor A, but a uniprot or refseq does exist,
            # that means the ecoli interactor A is a protein without orthologs, so continue to next interaction
            if (orthologs_A is None) & ((uniprot_A is not None) | (refseq_A is not None)): continue

            # same as for interactor A above
            if 'uniprotkb' in row['interactor_B']:
                uniprot_B = row['interactor_B'].split('uniprotkb:')[1].split('|')[0]
            if 'refseq' in row['interactor_B']:
                refseq_B = row['interactor_B'].split('refseq:')[1].split('|')[0]

            if uniprot_B is not None:
                orthologs_B = session.query(OrthologEcoli).filter_by(ortholog_uniprot = uniprot_B).all()
            if (orthologs_B is None) and (refseq_B is not None):
                orthologs_B = session.query(OrthologEcoli).filter_by(ortholog_refseq = refseq_B).all()
            if (orthologs_B is None) & ((uniprot_B is not None) | (refseq_B is not None)): continue

            # if both orthologs_A and orthologs_B are None, then there are no protein interactors for this
            # interaction, so move on to the next interaction
            if (orthologs_A is None) and (orthologs_B is None): continue

            # if there were no orthologs for interactor A (and no refseq or uniprot was found),
            # search the file for pubchem or chebi ids for interactor A (as it may be a metabolite)
            if orthologs_A is None:
                if 'chebi' in row['interactor_A']:
                    chebi = row['interactor_A'].split('CHEBI:')[1].split('|')[0][:-1]
                if 'pubchem' in row['altID_A']:
                    pubchem = row['altID_A'].split('pubchem:')[1].split('|')[0]
                if (chebi is None) & ('chebi' in row['altID_A']):
                    chebi = row['altID_A'].split('CHEBI:')[1].split('|')[0][:-1]
                # if no metabolite ids were found in the interaction row, then move on to the next interaction
                # because no interactor_A was identified
                if (chebi is None) & (pubchem is None): continue
                # if a pubchem or chebi id was found, then this interaction will be a p-m interaction, so
                # set the protein interactors(orthologs) to orthologs_B
                orthologs = orthologs_B
            # other case where orthologs_B were not identified so need to check if interactor B has metabolite ids
            elif orthologs_B is None:
                if 'chebi' in row['interactor_B']:
                    chebi = row['interactor_B'].split('CHEBI:')[1].split('|')[0][:-1]
                if 'pubchem' in row['altID_B']:
                    pubchem = row['altID_B'].split('pubchem:')[1].split('|')[0]
                if (chebi is None) & ('chebi' in row['altID_B']):
                    chebi = row['altID_B'].split('CHEBI:')[1].split('|')[0][:-1]
                if (chebi is None) & (pubchem is None): continue
                orthologs = orthologs_A

            # if one of the interactors was identified to be a metabolite, search for the metabolite and set metabolite
            # variable to that value. if the metabolite doesnt exist create it
            # Note: if this point was reached, it means one of the interactors had protein orthologs,
            # so we can safely create a new metabolite knowing it will have a protein interaction partner
            if (chebi is not None) | (pubchem is not None):
                id = None
                # preferentially set id for new metabolites to be chebi
                if chebi is not None:
                    id = chebi
                    metabolite = session.query(Metabolite).filter_by(chebi = chebi).first()
                # if no metabolite with chebi was found, but pubchem id exists, try to find
                # metabolite with that pubchem
                if (metabolite is None) & (pubchem is not None):
                    id = pubchem
                    metabolite = session.query(Metabolite).filter_by(pubchem = pubchem).first()
                # if no metabolite was found with pubchem or chebi id, create new metabolite
                if metabolite is None:
                    metabolite = Metabolite(id = id, chebi=chebi, pubchem=pubchem)
                    session.add(metabolite)
                # if a metabolite was found, update its chebi and pubchem if it has none
                else:
                    if metabolite.pubchem is None:
                        metabolite.pubchem = pubchem
                    if metabolite.chebi is None:
                        metabolite.chebi = chebi

            # list of interactor pairs for interaction
            interactors = []
            # if no metabolite was found for interaction, it is a p-p interaction, so iterate through
            # orthologs to create interactor pairs
            if metabolite is None:
                for ortholog_A in orthologs_A:
                    for ortholog_B in orthologs_B:
                        if (ortholog_A is not None) and (ortholog_B is not None):
                            # only add the interactor pair if the protein strains match
                            if ortholog_A.strain_protein == ortholog_B.strain_protein:
                                interactors.append([[ortholog_A.protein, ortholog_A.ortholog_id],
                                                    [ortholog_B.protein, ortholog_B.ortholog_id]])
            else:
                # if a metabolite was found, add pairs of all orthologs with metabolite to interactor pairs
                for ortholog in orthologs:
                    interactors.append([[metabolite, metabolite.id], [ortholog.protein, ortholog.ortholog_id]])

            # for each interactor pair, create interaction if it doesnt exist, otherwise update attributes
            for interactor_pair in interactors:
                homogenous = (interactor_pair[0][0] == interactor_pair[1][0])
                interaction = session.query(Interaction).filter(Interaction.interactors.contains(interactor_pair[0][0]),
                                                                Interaction.interactors.contains(interactor_pair[1][0]),
                                                                Interaction.homogenous == homogenous).first()
                if interaction is None:
                    # since one of the interactors may be a metabolite, set strain to match strain of protein
                    strain = None
                    if interactor_pair[0][0].type == 'p':
                        strain = interactor_pair[0][0].strain
                    else:
                        strain = interactor_pair[1][0].strain
                    # if interaction did not exist, set it to Ecoli ortholog derived
                    interaction = Interaction(strain=strain,
                                              interactors=[interactor_pair[0][0], interactor_pair[1][0]],
                                              type=(interactor_pair[0][0].type + '-' + interactor_pair[1][0].type),
                                              ortholog_derived='Ecoli')
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
                    ref_parameters = [comb[0][0], comb[0][1][0], comb[0][1][1], comb[0][1][2], comb[0][2],
                                      comb[1], comb[2]]
                    # dont add the ref_parametrs to the ref_parameter list if all parameters are None
                    if all(parameter is None for parameter in ref_parameters): continue
                    ref_parameter_list.append(ref_parameters)

                # in case the interaction already existed, make sure interactor_a and interactor_b variables for
                # new interaction reference match up with the first and second interactors of the existing
                # interaction (so it's easy to see which Pseudomonas interactor matches up with which Ecoli
                # ortholog)
                interactor_a, interactor_b = None, None
                if interaction.interactors[0] == interactor_pair[0][0]:
                    interactor_a = interactor_pair[0][1]
                    interactor_b = interactor_pair[1][1]
                else:
                    interactor_b = interactor_pair[0][1]
                    interactor_a = interactor_pair[1][1]

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
                        detection_method=ref[0], author_ln = ref[1], pub_date=ref[2], pmid=ref[3],
                        interaction_type=ref[4], source_db=ref[5], confidence=ref[6],
                        interactor_a=interactor_a, interactor_b=interactor_b).first()
                    # if nref doesn't exist, create and add it to the interaction's reference list,
                    # and add the source to the reference's sources
                    if nref is None:
                        nref = InteractionReference(
                            detection_method=ref[0], author_ln = ref[1], pub_date=ref[2], pmid=ref[3],
                            interaction_type=ref[4], source_db=ref[5], confidence=ref[6],
                            interactor_a=interactor_a, interactor_b=interactor_b)
                        interaction.references.append(nref)
                        nref.sources.append(nsource)
                    # if nref does exist, add the interaction and source to it's attributes if they aren't added
                    else:
                        if interaction not in nref.interactions:
                            nref.interactions.append(interaction)
                        if nsource not in nref.sources:
                            nref.sources.append(nsource)

    session.commit()
    print(source, session.query(Interaction).count())