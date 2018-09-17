if __name__ == '__main__':

    from Schema import Base
    from sqlalchemy import create_engine
    from sqlalchemy.orm import sessionmaker
    from Parsers.Parser import parse_all
    from DB_Statistics.DB_Statistics import get_db_stats_interactors, get_db_stats_interactions

    #('sqlite:///C:\\Users\\olgas\\Desktop\\PaIntDB.db')
    #'sqlite:////Users/olga/Desktop/PaIntDB.db'
    engine = create_engine('sqlite:////Users/olga/Desktop/PaIntDB.db')
    # if db is already created, next line can be commented out (tables do not have to be created again
    # to connect to the db)
    #Base.metadata.create_all(engine)
    Session = sessionmaker(bind=engine)
    session = Session()
    from Parsers.PAO1.STRING import parse
    #parse(session)
    # comment out whatever processes are unnecessary
    #parse_all(session)
    # get_db_stats_interactors('PAO1',session)
    #get_db_stats_interactions('PAO1',session)
    # get_db_stats_interactors('PA14', session)
    #get_db_stats_interactions('PA14', session)
    #from test import parse
    #parse()
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

    import csv
    from Schema import Interaction, InteractionSource, Interactor, InteractionReference, InteractionXref

    psimi_fields = ['ID(s) interactor A', 'ID(s) interactor B', 'Alt. ID(s) interactor A', 'Alt. ID(s) interactor B',
                    'Alias(es) interactor A', 'Alias(es) interactor B', 'Interaction detection method(s)',
                    'Publication 1st author(s)', 'Publication identifier(s)', 'Taxid interactor A',
                    'Taxid interactor B',
                    'Interaction type(s)', 'Source database(s)', 'Interaction identifier(s)', 'Confidence value(s)',
                    'Annotation(s) interactor A', 'Annotation(s) interactor B']

    filters = {'strain': ['PAO1'], 'interaction_type': ['p-p'], 'Ecoli_sources': ['None'], 'PAO1_sources':
        ['Geoff', 'XLinkDB', 'ADIPInteractomes(PAO1)', 'IMEx(PAO1)', 'IntAct(PAO1)',
                    'iRefIndex(PAO1)', 'mentha(PAO1)', 'MINT(PAO1)'], 'PA14_sources': ['None'],
               'verification': ['All']}
    sources_Ecoli = ['EcoCyc', 'RegulonDB(Ecoli)', 'IMEx(Ecoli)', 'BindingDB(Ecoli)', 'EBI-GOA-nonIntAct(Ecoli)',
                     'IntAct(Ecoli)', 'iRefIndex(Ecoli)', 'mentha(Ecoli)', 'MINT(Ecoli)', 'MPIDB(Ecoli)',
                     'UniProt(Ecoli)', 'DIP(Ecoli)', 'KEGG(Ecoli)']
    sources_PAO1 = ['Geoff', 'XLinkDB', 'Zhang', 'ADIPInteractomes(PAO1)', 'IMEx(PAO1)', 'IntAct(PAO1)',
                    'iRefIndex(PAO1)', 'mentha(PAO1)', 'MINT(PAO1)', 'Galan-Vasquez(PAO1)', 'KEGG(PAO1)']
    sources_PA14 = ['IMEx(PA14)', 'IntAct(PA14)', 'iRefIndex(PA14)', 'mentha(PA14)', 'MINT(PA14)', 'KEGG(PA14)',
                    'Galan-Vasquez(PA14)']
    filters_all = {'strain': ['PAO1', 'PA14'], 'interaction_type': ['p-p', 'p-m', 'p-bs'],
                   'Ecoli_sources': sources_Ecoli, 'PAO1_sources': sources_PAO1, 'PA14_sources': sources_PA14,
                   'verification': [0,1,2]}
    sources = []
    tfbs_sources = ['RegulonDB(Ecoli)', 'Galan-Vasquez(PAO1)']


    for filter in filters:
        if 'None' in filter:
            filters[filter] = ['None']
        elif 'All' in filter:
            filters[filter] = filters_all[filter]

    if filters['PAO1_sources'][0] != 'None':
        sources.append(filters['PAO1_sources'])
    if filters['PA14_sources'][0] != 'None':
        sources.append(filters['PA14_sources'])
    if filters['Ecoli_sources'][0] != 'None':
        sources.append(filters['Ecoli_sources'])
    if len(sources) == (len(filters_all['PAO1_sources']) + len(filters_all['PA14_sources']) +
                        len(filters_all['Ecoli_sources'])):
        sources = ['All']
    if len(sources) == 0:
        print("no sources selected")
    session = Session()
    print(sources)
    interactions = None
    if len(filters['strain']) ==2:
        interactions = session.query(Interaction)
    elif 'PAO1' in filters['strain']:
        interactions = session.query(Interaction).filter_by(strain='PAO1')
    elif 'PA14' in filters['strain']:
        interactions = session.query(Interaction).filter_by(strain='PA14')
    #
    # type = []
    #
    # if (len(filters['interaction_type']) == 1) & (filters['interaction_type'][0] == 'p-bs'):
    #     selected_tfbs = []
    #     if tfbs_sources[0] in sources:
    #         selected_tfbs.append(tfbs_sources[0])
    #     if tfbs_sources[1] in sources:
    #         selected_tfbs.append(tfbs_sources[1])
    #     if len(selected_tfbs) > 0:
    #         interactions = interactions.join(Interaction.sources).filter(InteractionSource.data_source.in_(selected_tfbs))
    # else:
    #     if ('p-p' in filters['interaction_type']) | ('p-bs' in filters['interaction_type']):
    #         type.append('p-p')
    #     if 'p-m' in filters['interaction_type']:
    #         type.append('p-m')
    #         type.append('m-p')
    #     if len(type) < 3:
    #         interactions = interactions.filter(Interaction.type.in_(type))
    #
    #     interactions = interactions.join(Interaction.sources)
    #     if sources[0] != 'All':
    #         interactions = interactions.filter(InteractionSource.data_source.in_(sources))
    #
    #
    # if len(filters['verification']) < 3:
    #     interactions = interactions.filter(InteractionSource.is_experimental.in_(filters['verification']))
    #
    #
    # interactions = interactions.join(Interaction.references)
    # interactions = interactions.join(Interaction.xrefs)

    file_writer = csv.DictWriter(open('output.csv', mode='x', newline=''), fieldnames=psimi_fields)
    file_writer.writeheader()
    for interaction in interactions.all():
        if interaction is None: continue
        interactor_ids, alt_ids, aliases = [], [], []
        is_protein = []

        for interactor in interaction.interactors:
            if interactor.name is None:
                aliases.append('')
            else:
                aliases.append(interactor.name)

            if interactor.type == 'p':
                is_protein.append(1)
                if interactor.uniprotkb == 'pc':
                    interactor_ids.append('uniprotkb:' + interactor.id)
                    alt_ids.append('')
                else:
                    interactor_ids.append('gene/locus_link:' + interactor.id)
                    alt_id = ''
                    # for xref in interactor.xrefs:
                    #     alt_id += xref.source + ':' + xref.accession + '|'
                    # if len(alt_id) != 0:
                    #     alt_id = alt_id[:-1]
                    alt_ids.append(alt_id)

            else:
                is_protein.append(0)
                id = None
                alt_id = ''
                if interactor.pubchem is not None:
                    id = 'pubchem:' + interactor.pubchem
                if interactor.chebi is not None:
                    chebi = 'chebi:"CHEBI:' + interactor.chebi + '"'
                    if id is None:
                        id = chebi
                    else:
                        alt_id += chebi + '|'
                if interactor.cas is not None:
                    cas = 'cas:' + interactor.cas
                    if id is None:
                        id = cas
                    else:
                        alt_id += cas + '|'
                if interactor.kegg is not None:
                    kegg = 'kegg:' + interactor.kegg
                    if id is None:
                        id = kegg
                    else:
                        alt_id += kegg + '|'
                if interactor.ecocyc is not None:
                    ecocyc = 'ecocyc:' + interactor.ecocyc
                    if id is None:
                        id = ecocyc
                    else:
                        alt_id += ecocyc + '|'

                if len(alt_id) != 0:
                    alt_id = alt_id[:-1]

                alt_ids.append(alt_id)

        if len(interactor_ids) == 1:
            interactor_ids.append(interactor_ids[0])
            alt_ids.append(alt_ids[0])
            aliases.append(aliases[0])
            is_protein.append(is_protein[0])

        taxid_A, taxid_B, taxid = None, None, None
        if interaction.strain == 'PAO1':
            taxid = 'taxid:208964(pseae)'
        else:
            taxid = 'taxid:208963(pseab)'
        if is_protein[0]:
            taxid_A = taxid
        if is_protein[1]:
            taxid_B = taxid

        # 0064 ortholog interaction (interologs mapping)
        refs = {'detection': [], 'author': [], 'pmid': [], 'type': [], 'db': [], 'xrefs': [],
                'confidence': [], 'annotations A': [], 'annotations B': []}
        author = []
        # sorted(list_with_none, key=lambda k: (k[col] is not None, k[col] != "", k[col]), reverse=True)
        for reference in interaction.references:
            refs['detection'].append(reference.detection_method)
            author.append([reference.author_ln, reference.pub_date])
            refs['type'].append(reference.interaction_type)
            refs['db'].append(reference.source_db)
            refs['confidence'].append(reference.confidence)
            refs['annotations A'].append(reference.interactor_a)
            refs['annotations B'].append(reference.interactor_b)

        for author_info in author:
            if author_info[0] is not None:
                author_info[0] += 'et al.'
                if author_info[1] is not None:
                    author_info[0] += ' (' + author_info[1] + ')'
            refs['author'].append(author_info[0])

        for xref in interaction.xrefs:
            if xref is None: continue
            refs['xrefs'].append(xref.data_source + ':' + xref.accession)
        for ref in refs:
            refs[ref] = ['-' if field is None else field for field in refs[ref]]
            if refs[ref].count(refs[ref][0]) == len(refs[ref]):
                refs[ref] = [refs[ref][0]]
            refs[ref] = ['|'.join(refs[ref])]

        file_writer.writerow({psimi_fields[0]: interactor_ids[0], psimi_fields[1]: interactor_ids[1],
                              psimi_fields[2]: alt_ids[0], psimi_fields[3]: alt_ids[1],
                              psimi_fields[4]: aliases[0], psimi_fields[5]: aliases[1],
                              psimi_fields[6]: refs['detection'][0], psimi_fields[7]: refs['author'][0],
                              psimi_fields[8]: refs['pmid'][0], psimi_fields[9]: taxid_A, psimi_fields[10]: taxid_B,
                              psimi_fields[11]: refs['type'][0], psimi_fields[12]: refs['db'][0],
                              psimi_fields[13]: refs['xrefs'][0], psimi_fields[14]: refs['confidence'][0],
                              psimi_fields[15]: refs['annotations A'], psimi_fields[15]: refs['annotations B']})