import csv
from DB_schema import Interactor, OrthologEcoli

# dicts for each strain of form:
# {strain 1 id: [strain 2 inparalogs]}
inparalogs_PAO1, inparalogs_PA14 = {}, {}

# function to go through process of mapping Ecoli and Pseudomonas orthologs from ortholuge file
def parse(session):
    # first parse inparalogs for both PAO1 and PA14 strains
    parse_inparalogs('Data/Ortholog/PAO1-Ecoli.csv', inparalogs_PAO1)
    parse_inparalogs('Data/Ortholog/PA14-Ecoli.csv', inparalogs_PA14)
    # after identifying inparalogs, parse orthologs for both strains
    parse_orthologs('Data/Ortholog/PAO1-Ecoli.csv', 'PAO1', inparalogs_PAO1, session)
    parse_orthologs('Data/Ortholog/PA14-Ecoli.csv', 'PA14', inparalogs_PA14, session)
    # fill in uniprot ids for ecoli orthologs
    parse_uniprot_ids('Data/Ecoli/Uniprot_IDs.csv', session)
    print('e_orthologs')


def parse_inparalogs(file, dict):
    # get all inparalogs and store in inparalogs dict
    with open(file) as csvfile:
        reader = csv.DictReader(csvfile)

        for row in reader:
            # ignore Non/Borderline SSD inparalogs (they will be ignored when parsing orthologs anyway)
            if (row['Ortholuge Class'] != '') & (row['Ortholuge Class'] != 'SSD'): continue

            #check for strain 1 (Ecoli) inparalogs
            if row['Strain 1 Inparalogs (Locus Tag/Name)'] != '':
                strain1_inparalogs = row['Strain 1 Inparalogs (Locus Tag/Name)'].split(']')

                # for all inparalogs, add them to given strain's inparalog dict
                for inparalog in strain1_inparalogs:
                    if inparalog == '': continue
                    inparalog_id = inparalog.split('[')[0]
                    if inparalog_id[0] == ';':
                        inparalog_id = inparalog_id[1:]
                    if inparalog_id in dict:
                        dict[inparalog_id].append(row['Locus Tag (Strain 2)'])
                    else:
                        dict[inparalog_id] = [row['Locus Tag (Strain 2)']]

            # check for strain 2 (Pseudomonas) inparalogs
            if row['Strain 2 Inparalogs (Locus Tag/Name)'] != '':
                strain2_inparalogs = row['Strain 2 Inparalogs (Locus Tag/Name)'].split(']')

                for inparalog in strain2_inparalogs:
                    if inparalog == '': continue
                    inparalog_id = inparalog.split('[')[0]
                    if inparalog_id[0] == ';':
                        inparalog_id = inparalog_id[1:]
                    if row['Locus Tag (Strain 1)'] in dict:
                        dict[row['Locus Tag (Strain 1)']].append(inparalog_id)
                    else:
                        dict[row['Locus Tag (Strain 1)']] = [inparalog_id]

# function to parse orthologs between Ecoli and Pseudomonas for given strain from ortholuge file,
# ignoring inparalogs
def parse_orthologs(file, strain, dict, session):
    # parse all orthologs (don't include inparalogs)
    with open(file) as csvfile:
        reader2 = csv.DictReader(csvfile)
        for row in reader2:
            # only care about RBBH (blank ortholuge class) or SSD orthologs
            if (row['Ortholuge Class'] != '') & (row['Ortholuge Class'] != 'SSD'): continue

            # don't include inparalog interactions
            if row['Locus Tag (Strain 1)'] in dict:
                if row['Locus Tag (Strain 2)'] in dict[row['Locus Tag (Strain 1)']]:
                    continue

            refseq = None
            if row['NCBI RefSeq Accession (Strain 1)'] != '':
                refseq = row['NCBI RefSeq Accession (Strain 1)']
            #check if pseudomonas protein exists in database, if it does then add the ecoli ortholog
            if session.query(Interactor).get(row['Locus Tag (Strain 2)']) is not None:
                ortholog = OrthologEcoli(protein_id=row['Locus Tag (Strain 2)'], strain_protein=strain,
                                         ortholog_id=row['Locus Tag (Strain 1)'], ortholog_refseq=refseq)

                if row['Ortholuge Class'] == '':
                    ortholog.ortholuge_classification = 'RBBH'
                else:
                    ortholog.ortholuge_classification = row['Ortholuge Class']
                session.add(ortholog)
    session.commit()

def parse_uniprot_ids(file, session):
    # get all uniprot ids and gene names for E. coli orthologs
    with open(file) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if row['Gene names'] != '':
                locus = row['Gene names'].split(' ')[0]
                name = None
                if row['Gene name'] != '':
                    name = row['Gene name']
                for ortholog in session.query(OrthologEcoli).filter_by(ortholog_id = locus).all():
                    if ortholog is not None:
                        ortholog.ortholog_uniprot = row['Entry']
                        ortholog.ortholog_name = row['Gene name']
    session.commit()