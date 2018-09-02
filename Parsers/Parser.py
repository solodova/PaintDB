if __name__ == '__main__':

    from Schema1 import Base, Interactor, Metabolite, Protein, InteractorXref, ProteinReference, OrthologPseudomonas, \
        OrthologEcoli, GeneOntology, Localization, InteractionReference, Interaction, InteractionSource
    from sqlalchemy import create_engine, or_
    from sqlalchemy.orm import sessionmaker
    import csv
    #/C:\\Users\\olgas\\Desktop\\PaIntDB.db
    engine = create_engine('sqlite://')
    Base.metadata.create_all(engine)
    Session = sessionmaker(bind=engine)
    session = Session()

    from Parsers.PAO1_PA14 import pseudomonas_db, Parse_PSIMI, regulatory_network
    from Parsers.PAO1 import Geoff_Winsor, STRING, xlinkdb, Zhang
    from Parsers.PAO1_PA14_Ecoli import KEGG
    from Parsers.Ecoli import EcoCyc, RegulonDB, ortholuge
    p1 = Protein(id='A')
    p2 = Protein(id='B')
    p3 = Protein(id='C')
    int1 = Interaction(type='1', interactors=[p1, p2], homogenous = 0)
    int2 = Interaction(type='2', interactors=[p1, p1], homogenous=1)
    int3 = Interaction(type='3', interactors=[p2, p3], homogenous=0)
    db1 = InteractionSource(data_source='db1')
    db2 = InteractionSource(data_source='db2')
    db3 = InteractionSource(data_source='db3')
    ref1 = InteractionReference(pmid='1')
    ref2 = InteractionReference(pmid='2')
    ref3 = InteractionReference(pmid='3')
    int1.references.append(ref1)
    int1.sources.append(db1)
    int1.sources.append(db2)
    int2.sources.append(db2)
    int3.sources.append(db3)
    session.add_all([int1, int2, int3, db1, db2, db3, p1, p2, p3])
    session.commit()
    #sources = ['db1', 'db2']
    query= session.query(Interaction).filter(Interaction.interactors.contains(p1),
                                             Interaction.interactors.contains(p2),
                                             Interaction.homogenous == 0).join(Interaction.sources, Interaction.references).first()

    print(query.sources, query.references)
    #query = session.query(Interaction).join(Interaction.sources).filter(InteractionSource.data_source.in_(sources)).filter(Interaction.type == '1')
    #query = session.query(Interaction).
    # for q in query.all():
    #     print(q.type)
    #     for s in q.sources:
    #         print(s.data_source)
    #     print(q.references)

    pseudomonas_db.parse_pseudomonasdb(session)

    #ortholuge.parse_ortholuge_ecoli(session)
    Geoff_Winsor.parse_geoff(session)
    print(session.query(Interaction).count())
    # xlinkdb.parse_xlinkdb(session)
    # Parse_PSIMI.parse_psimi_pseudomonas(session)
    # regulatory_network.parse_regulatory_network(session)
    # KEGG.get_kegg_compounds()
    # KEGG.parse_pseudomonas_kegg(session)
    # Zhang.parse_zhang(session)
    #ortholuge.parse_ortholuge(session)

    #print(session.query(Interaction).filter(Interaction.strain == 'PA14').count())
    #print(session.query(OrthologPseudomonas).filter(OrthologPseudomonas.strain_ortholog == 'PA14').count())
    #Parse_PSIMI.parse_ecoli_psimi(session)
    #EcoCyc.parse_ecocyc(session)
    #RegulonDB.parse_ecoli_regulondb(session)
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
