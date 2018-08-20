import csv

def parse():
    cols =[]
    all_num_occurrences = {}
    all_presence = {}
    with open('Data/Ecoli/PSICQUIC/BindingDB.txt') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        cols = reader.fieldnames


    for col in cols:
        with open('Data/Ecoli/PSICQUIC/BindingDB.txt') as csvfile:
            reader = csv.DictReader(csvfile, delimiter='\t')
            num_occurrences = {}

            for row in reader:
                if row[col] is None: continue
                fields = row[col].split('|')

                for id in fields:
                    if id.split(':')[0] in num_occurrences:
                        num_occurrences[id.split(':')[0]] += 1
                    else:
                        num_occurrences[id.split(':')[0]] = 1

            all_num_occurrences[col] = num_occurrences


    for col in cols:
        with open('Data/Ecoli/PSICQUIC/BindingDB.txt') as csvfile:
            reader = csv.DictReader(csvfile, delimiter='\t')
            presence = {}
            for row in reader:
                for id in all_num_occurrences[col]:
                    row_ids = []
                    for row_id in row[col].split('|'):
                        row_ids.append(row_id.split(':')[0])

                    if id in row_ids:
                        pass
                    else:
                        presence[id] = 'not present in all entries'
            all_presence[col] = presence

    for col in cols:
        print(col)
        print(all_num_occurrences[col])
        print(all_presence[col])
        print('')
