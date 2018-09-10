if __name__ == '__main__':

    from Schema1 import Base, Interactor, Metabolite, Protein, OrthologPseudomonas, \
        OrthologEcoli, GeneOntology, Localization, InteractionReference, Interaction, InteractionSource
    from sqlalchemy import create_engine, or_, not_, all_, any_, func
    from sqlalchemy.orm import sessionmaker
    import csv, datetime
    #('sqlite:///C:\\Users\\olgas\\Desktop\\PaIntDB.db')
    #'sqlite:////Users/olga/Desktop/PaIntDB.db'
    engine = create_engine('sqlite:///C:\\Users\\olgas\\Desktop\\PaIntDB.db')
    Base.metadata.create_all(engine)
    Session = sessionmaker(bind=engine)
    session = Session()

    from Parsers.PAO1_PA14 import pseudomonas_db, Parse_PSIMI_pseudomonas, regulatory_network, ortholuge_pseudomonas
    from Parsers.PAO1 import Geoff_Winsor, STRING, xlinkdb, Zhang
    from Parsers.PAO1_PA14_Ecoli import KEGG
    from Parsers.Ecoli import EcoCyc, RegulonDB, ortholuge_ecoli, Parse_PSIMI_ecoli
    from Schema2 import get_db_stats

    # i1 = Interaction(type = '1')
    # i2 = Interaction(type='2')
    # i3 = Interaction(type='3')
    # d1 = InteractionSource(data_source = 'a')
    # d2 = InteractionSource(data_source='b')
    # d3 = InteractionSource(data_source='c')
    # i1.sources.append(d1)
    # i1.sources.append(d2)
    # i2.sources.append(d1)
    # i3.sources.append(d2)
    # i3.sources.append(d3)
    # session.add_all([i1,i2,i3,d1,d2,d3])
    #
    # #d = session.query(InteractionSource)
    # # num PA14 = num with PA14 - num with PAO1 and PA14
    # # num ecoli = num with Ecoli - num with (PAO1 | PA14) and Ecoli
    # print(session.query(Interaction).join('sources').
    #                  group_by(Interaction).filter(Interaction.sources.any(InteractionSource.data_source.in_(['a']))).
    #                  filter(~(Interaction.sources.any(InteractionSource.data_source.in_(['b'])))).count())
    #filter(InteractionSource.data_source)
    #print(session.query(Protein).filter(Protein.strain == 'PAO1',Protein.uniprotkb == 'pc').count())
    #get_db_stats('PAO1')
    #get_db_stats('PA14')
    pseudomonas_db.parse(session)
    print(datetime.datetime.now())
    Geoff_Winsor.parse(session)
    print(datetime.datetime.now())
    xlinkdb.parse(session)
    print(datetime.datetime.now())
    Zhang.parse(session)
    print(datetime.datetime.now())
    regulatory_network.parse(session)
    print(datetime.datetime.now())
    Parse_PSIMI_pseudomonas.parse(session)
    print(datetime.datetime.now())
    KEGG.parse_pseudomonas(session)
    print(datetime.datetime.now())
    #print(session.query(Interaction).filter_by(strain = 'PAO1').count())
    #print(session.query(Interaction).filter_by(strain = 'PA14').count())
    ortholuge_pseudomonas.parse(session)
    print(datetime.datetime.now())

    #ortholuge_ecoli.parse(session)
    #print(datetime.datetime.now())
    KEGG.parse_ecoli(session)
    print(datetime.datetime.now())
    EcoCyc.parse(session)
    print(datetime.datetime.now())
    Parse_PSIMI_ecoli.parse(session)
    print(datetime.datetime.now())
    RegulonDB.parse(session)
    print(datetime.datetime.now())
    KEGG.update_metabolite_info_kegg(session)
    EcoCyc.update_metabolite_info_ecocyc(session)
    #print(session.query(Interaction).filter(Interaction.strain == 'PA14').count())
    #print(session.query(OrthologPseudomonas).filter(OrthologPseudomonas.strain_ortholog == 'PA14').count())
    # print(session.query(Interaction).filter(Interaction.strain == 'PAO1').count())
    # print(session.query(Interaction).filter(Interaction.strain == 'PAO1',
    #                                         Interaction.type == 'p-p').count())
    # print(session.query(Interaction).filter(Interaction.strain == 'PAO1',
    #                                         or_(Interaction.type == 'm-p',
    #                                             (Interaction.type == 'p-m'))).count())
    #
    # print(session.query(Interaction).filter(Interaction.strain == 'PA14').count())
    # print(session.query(Interaction).filter(Interaction.strain == 'PA14',
    #                                         Interaction.type == 'p-p').count())
    # print(session.query(Interaction).filter(Interaction.strain == 'PA14',
    #                                         or_(Interaction.type == 'm-p',
    #                                             (Interaction.type == 'p-m'))).count())
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
