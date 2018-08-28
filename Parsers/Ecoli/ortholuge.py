import csv
from Schema1 import Interactor, OrthologEcoli

inparalogs_PAO1 = {}
inparalogs_PA14 = {}

def parse_ortholuge_ecoli(session):
    parse_inparalogs('Orthologs/PAO1-Ecoli.csv', inparalogs_PAO1)
    parse_inparalogs('Orthologs/PAO1-Ecoli.csv', inparalogs_PA14)
    parse_orthologs('Orthologs/PAO1-Ecoli.csv', 'PAO1', inparalogs_PAO1, session)
    parse_orthologs('Orthologs/PA14-Ecoli.csv', 'PA14', inparalogs_PA14, session)
    parse_UniProt_IDs('Ecoli/Uniprot_IDs.csv', session)


def parse_inparalogs(file, dict):
    # get all inparalogs and store in inparalogs dict
    with open(file) as csvfile:
        reader = csv.DictReader(csvfile)

        for row in reader:
            # ignore Non/Borderline SSD
            if (row['Ortholuge Class'] != '') & (row['Ortholuge Class'] != 'SSD'): continue

            if row['Strain 1 Inparalogs (Locus Tag/Name)'] != '':
                strain1_inparalogs = row['Strain 1 Inparalogs (Locus Tag/Name)'].split(']')

                for inparalog in strain1_inparalogs:
                    if inparalog == '': continue
                    trimmed_inparalog = inparalog[:5]
                    if (inparalog[0] == ';'):
                        trimmed_inparalog = inparalog[1:6]
                    if (trimmed_inparalog in dict.keys()):
                        dict[trimmed_inparalog].append(row['Locus Tag (Strain 2)'])
                    else:
                        dict[trimmed_inparalog] = [row['Locus Tag (Strain 2)']]

            if (row['Strain 2 Inparalogs (Locus Tag/Name)'] != ''):
                strain2_inparalogs = row['Strain 2 Inparalogs (Locus Tag/Name)'].split(']')

                for inparalog in strain2_inparalogs:
                    if (inparalog == ''): continue
                    trimmed_inparalog = inparalog[:6]
                    if (inparalog[0] == ';'):
                        trimmed_inparalog = inparalog[1:7]
                    if (row['Locus Tag (Strain 1)'] in dict.keys()):
                        dict[row['Locus Tag (Strain 1)']].append(trimmed_inparalog)
                    else:
                        dict[row['Locus Tag (Strain 1)']] = [trimmed_inparalog]

def parse_orthologs(file, strain, dict, session):
    # parse all orthologs (don't include inparalogs)
    with open(file) as csvfile:
        reader2 = csv.DictReader(csvfile)
        for row in reader2:
            if (row['Ortholuge Class'] != '') & (row['Ortholuge Class'] != 'SSD'): continue

            if (row['Locus Tag (Strain 1)'] in dict.keys()):
                if (row['Locus Tag (Strain 2)'] in dict[row['Locus Tag (Strain 1)']]):
                    continue
            if (session.query(Interactor).filter(Interactor.id == row['Locus Tag (Strain 2)']).first() != None):
                ortholog = OrthologEcoli(protein_id=row['Locus Tag (Strain 2)'], strain_protein=strain,
                                         ortholog_id=row['Locus Tag (Strain 1)'], strain_ortholog='E. coli',
                                         ortholog_refseq=row['NCBI RefSeq Accession (Strain 1)'])

                if (row['Ortholuge Class'] == ''):
                    ortholog.ortholuge_classification = 'RBBH'
                else:
                    ortholog.ortholuge_classification = row['Ortholuge Class']
                session.add(ortholog)
        session.commit()

def parse_UniProt_IDs(file, session):
    # get all uniprot ids and gene names for E. coli orthologs
    with open(file) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if (row['Gene names'] != ''):
                locus = row['Gene names'].split(' ')[0]
                for ortholog in session.query(OrthologEcoli).filter(OrthologEcoli.ortholog_id == locus).all():
                    ortholog.ortholog_uniprot = row['Entry']
                    ortholog.ortholog_name = row['Gene name']
                session.commit()