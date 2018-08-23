if __name__ == '__main__':

    from Schema1 import Base, Interactor, Metabolite, Protein, InteractorXref, Reference, OrthologPseudomonas, \
        OrthologEcoli, GeneOntology, Localization
    from sqlalchemy import create_engine
    from sqlalchemy.orm import sessionmaker
    import csv

    engine = create_engine('sqlite:///:memory:', echo=True)
    Base.metadata.create_all(engine)
    Session = sessionmaker(bind=engine)
    session = Session()

    # stopped at chromatography tech
    experimental_psimi = ['MI:0045', 'MI:0401', 'MI:0400','MI:0004', 'MI:2288', 'MI:0019', 'MI:0006', 'MI:0007',
                          'MI:0402', 'MI:0225', 'MI:1218', 'MI:1028', 'MI:1029', 'MI:0858', 'MI:0946', 'MI:1017',
                          'MI:0729', 'MI:0096', 'MI:0963', 'MI:0676', 'MI:0008', 'MI:0225', 'MI:0081', 'MI:0089',
                          'MI:0678', 'MI:0695', 'MI:0092', 'MI:0946', 'MI:0095', 'MI:0921', 'MI:0405', 'MI:0034',
                          'MI:0009', 'MI:0073', 'MI:0084', 'MI:0048', 'MI:0899', 'MI:0900', 'MI:0066', 'MI:0108',
                          'MI:0098', 'MI:0115', 'MI:1087', 'MI:1031', 'MI:0813', 'MI:0440', 'MI:0892', 'MI:2189',
                          'MI:0947', 'MI:0411', 'MI:0047', 'MI:2283', 'MI:0049', 'MI:2167', 'MI:0657', 'MI:1232']
    def find_type(psi_code):
        return {
            '1110': 'predicted interaction',
            '2232': 'molecular association',
            '0914': 'association',
            '0915': 'physical association',
            '0407': 'direct interaction',
        }[psi_code]

    def is_experimental_psimi(psi_code):
        return psi_code in experimental_psimi


    from File_Descriptions.file_desc import parse
    parse()

    # info = {}
    # info_num = {"Locus Tag": 0, "Name": 0, "Product Name": 0, "NCBI Accession": 0}
    # with open('Data/PAO1/Gene_Info/pseudomonas_db_info.csv') as csvfile:
    #     fieldnames = ["Sequence", "Locus Tag", "Feature Type", "Start", "End", "Strand", "Name", "Product Name",
    #                   "Synonyms", "NCBI Accession"]
    #     interest_fields = ["Locus Tag", "Name", "Product Name", "NCBI Accession"]
    #     reader = csv.DictReader(csvfile, fieldnames=fieldnames)
    #     row_num = 0
    #     for row in reader:
    #         if (row['Feature Type'] == 'CDS') & (row_num >= 3):
    #             info[row['Locus Tag'][:-1]] = 1
    #             for field in interest_fields:
    #                 if row[field] != '':
    #                     info_num[field] += 1
    #         row_num += 1
    #     print(info_num)
    #
    # xrefs = {}
    # num_xrefs = {}
    # num_xrefs_unique = {}
    # with open('Data/PAO1/Gene_Info/pseudomonas_db_xrefs.csv') as csvfile:
    #     types = ['RefSeq Accession', 'UniProtKB Accession', 'UniProtKB ID', 'GI Number', 'Uniparc', 'UniRef100 ID',
    #              'UniRef90 ID', 'UniRef50 ID']
    #
    #     #types = ['UniProtKB Accession', 'RefSeq Accession']
    #     for locus in info:
    #         xrefs[locus] = {}
    #         for type in types:
    #             xrefs[locus][type] = 0
    #             num_xrefs[type] = 0
    #
    #     reader = csv.DictReader(csvfile)
    #     for row in reader:
    #         if (row['Locus Tag'] in info.keys()):
    #             for type in types:
    #                 if row[type] != '':
    #                     num_xrefs[type] += 1
    #                     xrefs[row['Locus Tag']][type] += 1
    #
    #     for type in types:
    #         num_xrefs_unique[type] = 0
    #
    #     for locus in xrefs:
    #         for type in xrefs[locus]:
    #             if xrefs[locus][type] >= 1:
    #                 num_xrefs_unique[type] += 1
    #
    #     print(num_xrefs)
    #     print(num_xrefs_unique)
    #
    #
    # localizations = {}
    # num_localizations = 0
    # num_localizations_unique = 0
    # for locus in info:
    #     localizations[locus] = 0
    # with open('Data/PAO1/Gene_Info/pseudomonas_db_localizations.csv') as csvfile:
    #     reader = csv.DictReader(csvfile)
    #     for row in reader:
    #         if row['Locus Tag'] in info:
    #             localizations[row['Locus Tag']] += 1
    #             num_localizations += 1
    #
    #     for locus in localizations:
    #         if localizations[locus] >= 1:
    #             num_localizations_unique += 1
    #
    #     print(num_localizations, num_localizations_unique)
    #
    # ontologies = {}
    # num_ontologies = 0
    # num_ontologies_unique = 0
    # for locus in info:
    #     ontologies[locus] = 0
    # with open('Data/PAO1/Gene_Info/pseudomonas_db_go.csv') as csvfile:
    #     reader = csv.DictReader(csvfile)
    #     for row in reader:
    #         if row['Locus Tag'] in info:
    #             ontologies[row['Locus Tag']] += 1
    #             num_ontologies += 1
    #
    #     for locus in ontologies:
    #         if ontologies[locus] >= 1:
    #             num_ontologies_unique += 1
    #
    #     print(num_ontologies, num_ontologies_unique)
    #
