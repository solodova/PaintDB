import csv, os


input_prefix = 'Data/'
output_prefix = 'File_Descriptions/'

def parse():
    cols = []
    all_num_occurrences = {}
    all_presence = {}
    sample_entry = {}
    num_elems = {}
    file_names = ['Ecoli/PSICQUIC/BindingDB.txt', 'Ecoli/PSICQUIC/ChEMBL.txt','Ecoli/PSICQUIC/EBI-GOA-nonIntAct.txt',
                  'Ecoli/PSICQUIC/IMEx.txt', 'Ecoli/PSICQUIC/IntAct.txt', 'Ecoli/PSICQUIC/iRefIndex.txt',
                  'Ecoli/PSICQUIC/mentha.txt', 'Ecoli/PSICQUIC/MINT.txt', 'Ecoli/PSICQUIC/MPIDB.txt',
                  'Ecoli/PSICQUIC/UniProt.txt', 'Ecoli/DIP.txt', 'Ecoli/IntAct.txt',
                  'PA14/PSICQUIC/IMEx.txt', 'PA14/PSICQUIC/iRefIndex.txt', 'PA14/PSICQUIC/mentha.txt',
                  'PA14/PSICQUIC/MINT.txt', 'PA14/IntAct.txt',
                  'PAO1/PSICQUIC/ADIPInteractomes.txt', 'PAO1/PSICQUIC/IMEx.txt', 'PAO1/PSICQUIC/IntAct.txt',
                  'PAO1/PSICQUIC/iRefIndex.txt', 'PAO1/PSICQUIC/mentha.txt', 'PAO1/PSICQUIC/MINT.txt',
                  'PAO1/PSICQUIC/MPIDB.txt', 'PAO1/IntAct.txt']


    #for file in file_names:
    file = 'Ecoli/PSICQUIC/UniProt.txt'

    input_file = input_prefix+file
    with open(input_file) as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        cols = reader.fieldnames

    for col in cols:
        num_elems[col] = {}
        all_num_occurrences[col] = {}
        with open(input_file) as csvfile:
            reader = csv.DictReader(csvfile, delimiter='\t')
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
            reader= csv.DictReader(csvfile, delimiter='\t')
            for row in reader:
                if row[col] is None: continue
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
            reader = csv.DictReader(csvfile, delimiter='\t')
            presence_prefix = 'not present in all entries: '
            presence = {}
            for row in reader:
                for id in all_num_occurrences[col]:
                    row_ids = []
                    for row_id in row[col].split('|'):
                        row_ids.append(row_id.split(':')[0])

                    if id in row_ids:
                        pass
                    else:
                        if id not in presence:
                            presence[id] = 1

            all_presence[col] = presence_prefix+(', '.join(presence.keys()))

    output_file = output_prefix+file
    new_file = open(output_file, mode='x')
    for col in cols:
        new_file.write('column name: '+col+'\n')
        if col in sample_entry:
            new_file.write('sample entries: ' + str(sample_entry[col]) + '\n')
        else:
            new_file.write('sample entries: None\n')
        new_file.write(str(all_num_occurrences[col])+'\n')
        new_file.write(all_presence[col] + '\n')
        new_file.write('num elements in column: ' + str(num_elems[col]) + '\n\n')

    new_file.close()

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