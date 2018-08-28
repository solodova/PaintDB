import csv, os


input_prefix = 'Data/'
output_prefix = 'File_Descriptions/'

def parse():
    # dictionary of number of times each type of id is present for each column of given file
    num_occurences_ids = {}
    #dictionary of lists of ids not present in each entry for each column of given file
    ids_not_present_in_all = {}
    #dictionary of number of times given element can appear in one folumn entry
    num_appearances_in_entry = {}
    #dictionary of sample entries for each column of given file
    sample_entry = {}
    # dictionary of number of separate fields present in each column of given file
    num_fields_in_col = {}
    file_names = ['Ecoli/PSICQUIC/BindingDB.txt', 'Ecoli/PSICQUIC/ChEMBL.txt','Ecoli/PSICQUIC/EBI-GOA-nonIntAct.txt',
                  'Ecoli/PSICQUIC/IMEx.txt', 'Ecoli/PSICQUIC/IntAct.txt', 'Ecoli/PSICQUIC/iRefIndex.txt',
                  'Ecoli/PSICQUIC/mentha.txt', 'Ecoli/PSICQUIC/MINT.txt', 'Ecoli/PSICQUIC/MPIDB.txt',
                  'Ecoli/PSICQUIC/UniProt.txt', 'Ecoli/DIP.txt', 'Ecoli/IntAct.txt',
                  'PA14/PSICQUIC/IMEx.txt', 'PA14/PSICQUIC/iRefIndex.txt', 'PA14/PSICQUIC/mentha.txt',
                  'PA14/PSICQUIC/MINT.txt', 'PA14/IntAct.txt',
                  'PAO1/PSICQUIC/ADIPInteractomes.txt', 'PAO1/PSICQUIC/IMEx.txt', 'PAO1/PSICQUIC/IntAct.txt',
                  'PAO1/PSICQUIC/iRefIndex.txt', 'PAO1/PSICQUIC/mentha.txt', 'PAO1/PSICQUIC/MINT.txt',
                  'PAO1/PSICQUIC/MPIDB.txt', 'PAO1/IntAct.txt']
    cols = ['interactor_A', 'interactor_B', 'altID_A', 'altID_B', 'alias_A', 'alias_B', 'detection',
            'publication', 'publication_ID', 'taxid_A', 'taxid_B', 'type', 'source_db', 'identifier',
            'confidence']
    of_interest_ecoli_cols = ['interactor_A', 'interactor_B', 'altID_A', 'altID_B', 'detection',
            'publication', 'publication_ID', 'type', 'source_db', 'identifier', 'confidence']
    #for file in file_names:
    file = file_names[24]

    input_file = input_prefix+file

    for col in of_interest_ecoli_cols:
        num_fields_in_col[col] = {}
        num_occurences_ids[col] = {}
        with open(input_file) as csvfile:
            reader = csv.DictReader(csvfile, delimiter='\t', fieldnames=cols)
            next(reader)
            for row in reader:
                fields = []
                if row[col] is None:
                    fields.append('None')
                elif row[col] == '':
                    fields.append('empty str')
                else:
                    fields = row[col].split('|')

                num_fields = 0
                for field in fields:
                    if field != '':
                        num_fields += 1

                if num_fields not in num_fields_in_col[col]:
                    num_fields_in_col[col][num_fields] = 1
                else:
                    num_fields_in_col[col][num_fields] += 1

                for field in fields:
                    id = field.split(':')[0]
                    if ((col != 'publication') & (':' in field)) | \
                            (id == '-') | (id == 'None') | (id == 'empty str') :
                        if id in num_occurences_ids[col]:
                            num_occurences_ids[col][id] += 1
                        else:
                            num_occurences_ids[col][id] = 1

    for col in of_interest_ecoli_cols:
        sample_entry[col] = []
        col_ids_todo = list(num_occurences_ids[col].keys())
        publication = []
        with open(input_file) as csvfile:
            reader= csv.DictReader(csvfile, delimiter='\t', fieldnames=cols)
            next(reader)
            for row in reader:
                fields = []
                if row[col] is None:
                    fields.append('None')
                elif row[col] == '':
                    fields.append('empty str')
                else:
                    fields = row[col].split('|')

                for field in fields:
                    if len(field.split(':')) == 1:
                        if (field == '-') | (field == 'None') | (field == 'empty str'):
                            if field in col_ids_todo:
                                sample_entry[col].append(field)
                                col_ids_todo.remove(field)
                        else:
                            if col == 'publication':
                                if len(publication) == 0:
                                    publication.append(field)
                                    sample_entry[col].append(field)
                                elif len(publication) == 1:
                                    if (('(' in publication[0]) and ('(' not in field)) | \
                                        (('(' not in publication[0]) and ('(' in field)):
                                        publication.append(field)
                                        sample_entry[col].append(field)
                    else:
                        if field.split(':')[0] in col_ids_todo:
                            sample_entry[col].append(field)
                            col_ids_todo.remove(field.split(':')[0])

    for col in of_interest_ecoli_cols:
        num_appearances_in_entry[col] ={}
        ids_not_present_in_all[col] = {}
        with open(input_file) as csvfile:
            reader = csv.DictReader(csvfile, delimiter='\t', fieldnames=cols)
            next(reader)
            for row in reader:
                row_ids = []
                if row[col] is None:
                    row_ids.append('None')
                elif row[col] == '':
                    row_ids.append('empty str')
                else:
                    for row_id in row[col].split('|'):
                        row_ids.append(row_id.split(':')[0])

                for id in num_occurences_ids[col]:
                    if (id not in row_ids) & \
                            ((col != 'publication') | (id == '-') | (id == 'None') | (id == 'empty str')):
                        if id not in ids_not_present_in_all[col]:
                            ids_not_present_in_all[col][id] = 1
                        else:
                            ids_not_present_in_all[col][id] += 1

                for id_option in num_occurences_ids[col]:
                    row_matching_ids = 0
                    for row_id in row_ids:
                        if row_id == id_option:
                            row_matching_ids += 1

                    if id_option in num_appearances_in_entry[col]:
                        if row_matching_ids in num_appearances_in_entry[col][id_option]:
                            num_appearances_in_entry[col][id_option][row_matching_ids] += 1
                        else:
                            num_appearances_in_entry[col][id_option][row_matching_ids] = 1
                    else:
                        num_appearances_in_entry[col][id_option] = {}
                        num_appearances_in_entry[col][id_option][row_matching_ids] = 1

    output_file = output_prefix+file
    new_file = open(output_file, mode='x')
    for col in of_interest_ecoli_cols:
        new_file.write('column name: '+col+'\n')
        new_file.write('sample entries: ' + str(sample_entry[col]) + '\n')
        new_file.write('total num occurences of each id:' + str(num_occurences_ids[col])+'\n')
        new_file.write( 'not present in all entries: ' + str(ids_not_present_in_all[col]) + '\n')
        new_file.write('number of each id in entries:' + str(num_appearances_in_entry[col]) + '\n')
        new_file.write('num elements in column: ' + str(num_fields_in_col[col]) + '\n\n')

    new_file.close()
    ############

    # #for file in file_names:
    # file = 'Ecoli/PSICQUIC/UniProt.txt'
    #
    # input_file = input_prefix+file
    # with open(input_file) as csvfile:
    #     reader = csv.DictReader(csvfile, delimiter='\t', fieldnames=cols)
    #
    # for col in cols:
    #     num_elems[col] = {}
    #     all_num_occurrences[col] = {}
    #     with open(input_file) as csvfile:
    #         reader = csv.DictReader(csvfile, delimiter='\t')
    #         for row in reader:
    #
    #             if row[col] is None: continue
    #             fields = row[col].split('|')
    #
    #             num_elem = 0
    #             for field in fields:
    #                 if field != '':
    #                     num_elem += 1
    #
    #             if num_elem not in num_elems[col]:
    #                 num_elems[col][num_elem] = 1
    #             else:
    #                 num_elems[col][num_elem] += 1
    #
    #             for id in fields:
    #                 if len(id.split(':')) == 1:
    #                     if id == '-':
    #                         if id in all_num_occurrences[col]:
    #                             all_num_occurrences[col][id] += 1
    #                         else:
    #                             all_num_occurrences[col][id] = 1
    #                 else:
    #                     if id.split(':')[0] in all_num_occurrences[col]:
    #                         all_num_occurrences[col][id.split(':')[0]] += 1
    #                     else:
    #                         all_num_occurrences[col][id.split(':')[0]] = 1
    #
    # for col in cols:
    #     sample_entry[col] = []
    #     ids = list(all_num_occurrences[col].keys())
    #     with open(input_file) as csvfile:
    #         reader= csv.DictReader(csvfile, delimiter='\t')
    #         for row in reader:
    #             if row[col] is None: continue
    #             fields = row[col].split('|')
    #
    #             for id in fields:
    #                 if len(id.split(':')) == 1:
    #                     if id == '-':
    #                         if id in ids:
    #                             sample_entry[col].append(id)
    #                             ids.remove(id)
    #                     else:
    #                         if (len(all_num_occurrences[col].keys()) + 1) > len(sample_entry[col]):
    #                             sample_entry[col].append(id)
    #                 else:
    #                     if id.split(':')[0] in ids:
    #                         sample_entry[col].append(id)
    #                         ids.remove(id.split(':')[0])
    #
    # for col in cols:
    #     with open(input_file) as csvfile:
    #         reader = csv.DictReader(csvfile, delimiter='\t')
    #         presence_prefix = 'not present in all entries: '
    #         presence = {}
    #         for row in reader:
    #             for id in all_num_occurrences[col]:
    #                 row_ids = []
    #                 for row_id in row[col].split('|'):
    #                     row_ids.append(row_id.split(':')[0])
    #
    #                 if id in row_ids:
    #                     pass
    #                 else:
    #                     if id not in presence:
    #                         presence[id] = 1
    #
    #         all_presence[col] = presence_prefix+(', '.join(presence.keys()))
    #
    # output_file = output_prefix+file
    # new_file = open(output_file, mode='x')
    # for col in cols:
    #     new_file.write('column name: '+col+'\n')
    #     if col in sample_entry:
    #         new_file.write('sample entries: ' + str(sample_entry[col]) + '\n')
    #     else:
    #         new_file.write('sample entries: None\n')
    #     new_file.write(str(all_num_occurrences[col])+'\n')
    #     new_file.write(all_presence[col] + '\n')
    #     new_file.write('num elements in column: ' + str(num_elems[col]) + '\n\n')
    #
    # new_file.close()

def parse_string():
    cols = []
    all_num_occurrences = {}
    all_presence = {}
    sample_entry = {}
    num_elems = {}
    file = 'PAO1/STRING.txt'
    cols = ['interactor_A', 'interactor_B', 'altID_A', 'altID_B', 'alias_A', 'alias_B', 'detection',
            'publication', 'publication_ID', 'taxid_A', 'taxid_B', 'type', 'source_db', 'identifier',
            'confidence']

    input_file = input_prefix + file
    with open(input_file) as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t', fieldnames=cols)

    for col in cols:
        num_elems[col] = {}
        all_num_occurrences[col] = {}
        with open(input_file) as csvfile:
            reader = csv.DictReader(csvfile, delimiter='\t', fieldnames=cols)
            for row in reader:

                if row[col] is None: continue
                fields = row[col].split('|')

                num_elem = 0
                for field in fields:
                    if field != '':
                        num_elem += 1

                if num_elem not in num_elems[col]:
                    num_elems[col][num_elem] = 1
                else:
                    num_elems[col][num_elem] += 1

                for id in fields:
                    if len(id.split(':')) == 1:
                        if id == '-':
                            if id in all_num_occurrences[col]:
                                all_num_occurrences[col][id] += 1
                            else:
                                all_num_occurrences[col][id] = 1
                    else:
                        if id.split(':')[0] in all_num_occurrences[col]:
                            all_num_occurrences[col][id.split(':')[0]] += 1
                        else:
                            all_num_occurrences[col][id.split(':')[0]] = 1

    for col in cols:
        sample_entry[col] = []
        ids = list(all_num_occurrences[col].keys())
        with open(input_file) as csvfile:
            reader = csv.DictReader(csvfile, delimiter='\t', fieldnames=cols)
            for row in reader:
                fields = row[col].split('|')

                for id in fields:
                    if len(id.split(':')) == 1:
                        if id == '-':
                            if id in ids:
                                sample_entry[col].append(id)
                                ids.remove(id)
                        else:
                            if (len(all_num_occurrences[col].keys()) + 1) > len(sample_entry[col]):
                                sample_entry[col].append(id)
                    else:
                        if id.split(':')[0] in ids:
                            sample_entry[col].append(id)
                            ids.remove(id.split(':')[0])

    for col in cols:
        with open(input_file) as csvfile:
            reader = csv.DictReader(csvfile, delimiter='\t', fieldnames=cols)
            presence_prefix = 'not present in all entries: '
            presence = {}
            for row in reader:
                for id in all_num_occurrences[col]:
                    row_ids = []
                    for row_id in row[col].split('|'):
                        row_ids.append(row_id.split(':')[0])

                    if id not in row_ids:
                        if id not in presence:
                            presence[id] = 1

            all_presence[col] = presence_prefix + (', '.join(presence.keys()))

    output_file = output_prefix + file
    new_file = open(output_file, mode='x')
    for col in cols:
        new_file.write('column name: ' + col + '\n')
        if col in sample_entry:
            new_file.write('sample entries: ' + str(sample_entry[col]) + '\n')
        else:
            new_file.write('sample entries: None\n')
        new_file.write(str(all_num_occurrences[col]) + '\n')
        new_file.write(all_presence[col] + '\n')
        new_file.write('num elements in column: ' + str(num_elems[col]) + '\n\n')

    new_file.close()