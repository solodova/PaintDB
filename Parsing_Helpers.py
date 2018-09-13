import itertools

experimental_psimi = ['0045', '0401', '0013', '0254', '0428', '1088', '0255', '0090', '0400', '1232', '0091', '0027',
                      '0030', '0982', '0415', '0417', '1036', '0953', '2197', '0968', '0016', '0943', '1311', '0894',
                      '0043', '1219', '1086', '0928', '0051', '0964', '0859', '0065', '0067', '2169', '1247', '0071',
                      '0893', '0891', '0077', '0938', '0099', '0888', '1235', '0966', '0114', '2224', '0439', '0441',
                      '0872', '0663', '0040', '0416', '0426', '2213', '0827', '1089', '1192', '0257', '2285', '0256',
                      '0014', '0010', '0011', '0809', '0111', '0231', '0895', '0097', '1203', '0370', '0232', '0004',
                      '0008', '0405', '0034', '1087', '1031', '0813', '0440', '0892', '0657', '0226', '0227', '0028',
                      '0029', '1022', '1211', '0430', '0031', '0807', '0983', '0889', '1005', '0989', '1142', '1147',
                      '1137', '0990', '1309', '0406', '2281', '0984', '2216', '1138', '0996', '1006', '0870', '1011',
                      '1009', '0998', '1026', '0999', '1007', '1000', '1249', '1001', '0515', '1010', '1034', '0879',
                      '0979', '0434', '1145', '0972', '0841', '0696', '1325', '1008', '0997', '1229', '0602', '0605',
                      '1191', '0949', '2198', '0969', '1342', '2196', '1038', '0107', '0069', '0944', '0041', '0042',
                      '0905', '0229', '0012', '0017', '1030', '0052', '0053', '1016', '0054', '0055', '0510', '0976',
                      '0965', '0038', '0104', '0826', '2170', '2171', '1104', '1103', '0425', '0825', '0824', '0410',
                      '1024', '0020', '1037', '1204', '0655', '0369', '1320', '0432', '0726', '0588', '0018', '2288',
                      '0019', '0096', '0676', '0225', '0081', '0089', '0921', '0073', '0084', '0098', '0115', '2189',
                      '0947', '0411', '0047', '0049', '2167', '0729', '2191', '2194', '1312', '0404', '0808', '0413',
                      '0887', '1354', '0991', '0435', '0508', '1252', '1236', '1003', '1002', '1004', '0516', '1035',
                      '0920', '0880', '0419', '0514', '1019', '0424', '0697', '0698', '0699', '0700', '0603', '1190',
                      '0604', '0901', '1189', '1183', '0814', '1352', '1313', '2199', '1246', '1238', '0009', '0420',
                      '0509', '0511', '1321', '0112', '0437', '0438', '0728', '0727', '0916', '0397', '0398', '1356',
                      '0006', '0007', '0402', '0858', '0946', '1017', '0963', '0678', '0092', '0095', '0048', '0066',
                      '0108', '2283', '2192', '2193', '2188', '2195', '0276', '0412', '0992', '0993', '0994', '0995',
                      '0513', '0512', '0423', '0606', '2168', '1184', '1314', '1113', '1111', '1218', '1028', '1029',
                      '0695', '0899', '0900', '1187', '2277', '2215', '1112', '0399']


# function to test whether psimi detection code is experimental
def is_experimental_psimi(psi_code):
    return psi_code in experimental_psimi

#function which returns whether or not there is an experimental psimi detection in the given interaction row
# Note: in order to work, cols for the row being passed in must be specfied to defaults in Parse_PSIMI.py files
def is_experimental_interaction(row):
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

    return is_experimental


# function to obtain list of all reference parameters present for a given interaction row (psimi formatted)
# Note: in order to work, cols for the row must be specfied to defaults in Parse_PSIMI.py files
def get_psimi_ref_list(row):
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
        # remove trailing comma from author name
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

    return ref_parameter_list