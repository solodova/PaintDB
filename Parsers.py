if __name__ == '__main__':

    import csv
    from Bio.KEGG.REST import kegg_list, kegg_get, kegg_conv
    from Bio.KEGG.KGML.KGML_parser import read
    from Schema3 import Interactor, Interaction, InteractionReference, InteractorXref, Reference, \
        Localization, GeneOntology, InteractionXref, OrthologPseudomonas, OrthologEcoli, Base
    from sqlalchemy import create_engine, or_
    from sqlalchemy.orm import sessionmaker

    engine = create_engine('sqlite:///:memory:', echo=True)
    Base.metadata.create_all(engine)
    Session = sessionmaker(bind=engine)

    session = Session()


    def find_type(psi_code):
        return {
            '1110': 'predicted interaction',
            '2232': 'molecular association',
            '0914': 'association',
            '0915': 'physical association',
            '0407': 'direct interaction',
        }[psi_code]


    pathways = kegg_list(database='pathway', org='pae').read().split('path:')
    pathwaylist = []
    for path in pathways:
        if (path != ''):
            pathwaylist.append(path[:8])
    # print(len(pathwaylist))

    jsonpaths = []
    for path in pathwaylist:
        jsonpaths.append(read(kegg_get(option='kgml', dbentries=path)))

    print(len(jsonpaths))
    jsonpaths[0].read()


    def getGeneInfo(file, strain):
        with open(file) as csvfile:
            fieldnames = ["Sequence", "Locus Tag", "Feature Type", "Start", "End", "Strand", "Name", "Product Name",
                          "Synonyms", "Accession", "GI", "Length (nucleotides)", "Isoelectric Point", "MW",
                          "Hydropathicity", "Charge at pH 7", "Length (amino acids)", "Nucleotide Sequence",
                          "Amino Acid Sequence"]
            reader = csv.DictReader(csvfile, fieldnames=fieldnames)
            row_num = 0
            string = []
            for row in reader:
                if (row['Feature Type'] == 'CDS') & (row_num >= 3):
                    interactor = Interactor(id=row['Locus Tag'][:-1], name=row['Name'], type=row['Feature Type'],
                                            strain=strain, description=row['Product Name'], accession=row['Accession'],
                                            start_site=row['Start'], end_site=row['End'])
                    session.add(interactor)
                    string.append(row['Locus Tag'][:-1] + ':' + row['Name'])
                row_num += 1
            session.commit()


    def getXrefs(file, strain):
        with open(file) as csvfile:
            types = ['RefSeq Accession', 'UniProtKB Accession', 'UniProtKB ID', 'GI Number', 'Uniparc', 'UniRef100 ID',
                     'UniRef90 ID', 'UniRef50 ID']
            reader = csv.DictReader(csvfile)
            for row in reader:
                if (row['Feature Type'] == 'CDS'):
                    for type in types:
                        if (row[type] != ''):
                            if (session.query(InteractorXref).filter((InteractorXref.accession == row[type]),
                                                                     (InteractorXref.interactor_id == row[
                                                                         'Locus Tag'])).count() == 0):
                                xref = InteractorXref(accession=row[type], interactor_id=row['Locus Tag'],
                                                      data_source=type)
                                session.add(xref)

                            if (type == 'UniProtKB Accession'):
                                if (session.query(Interactor).filter(
                                        Interactor.id == row['Locus Tag']).first() != None):
                                    interactor = session.query(Interactor).filter(
                                        Interactor.id == row['Locus Tag']).one()
                                    interactor.uniprotkb = row[type]
                                    session.commit()
                                complex_interactors = \
                                    session.query(Interactor).filter(
                                        Interactor.uniprotkb == row[type]).count()
                                if (complex_interactors == 2):
                                    interactor = Interactor(id=row[type], strain=strain,
                                                            type='protein complex')
                                    session.add(interactor)
                                    session.commit()
            session.commit()


    def getLocalizations(file):
        with open(file) as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                localization = Localization(interactor_id=row['Locus Tag'],
                                            localization=row['Subcellular Localization'],
                                            confidence=row['Confidence'])
                session.add(localization)
            session.commit()


    def getOntology(file):
        with open(file) as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                # if (session.query(Interactor).filter(Interactor.id == row['Locus Tag']).first() != None):
                ontology = GeneOntology(interactor_id=row['Locus Tag'], accession=row['Accession'],
                                        go_term=row['GO Term'],
                                        ontology=row['Namespace'], evidence_code=row['GO Evidence Code'],
                                        pmid=row['PMID'],
                                        eco_code=row['Evidence Ontology ECO Code'],
                                        eco_term=row['Evidence Ontology Term'])
                session.add(ontology)
            session.commit()


    getGeneInfo('Pseudomonas_aeruginosa_PAO1_107.csv', 'PAO1')
    getGeneInfo('Pseudomonas_aeruginosa_UCBPP-PA14_109.csv', 'PA14')

    getXrefs('pseudomonas_db_PAO1.csv', 'PAO1')
    getXrefs('pseudomonas_db_PA14.csv', 'PA14')

    getLocalizations('localizations_in_Pseudomonas_aeruginosa_PAO1_(Reference).csv')
    getLocalizations('localizations_in_Pseudomonas_aeruginosa_UCBPP-PA14.csv')

    getOntology('gene_ontology_csv_PAO1.csv')
    getOntology('gene_ontology_csv_PA14.csv')


    def parse_IMEx(file, strain, taxid):
        with open(file) as csvfile:
            reader = csv.DictReader(csvfile, delimiter='\t')
            for row in reader:
                interactors = []

                if ((row['Taxid interactor A'].split('|')[0] == taxid) &
                        (row['Taxid interactor B'].split('|')[0] == taxid)):

                    if (session.query(Interactor).filter(
                            Interactor.id == row['#ID(s) interactor A'].split(':')[1]).first() != None):
                        interactors.append(session.query(Interactor).filter(
                            Interactor.id == row['#ID(s) interactor A'].split(':')[1]).one())
                    elif (session.query(Interactor).filter(
                                Interactor.uniprotkb == row['#ID(s) interactor A'].split(':')[1]).first() != None):
                            interactors.append(session.query(Interactor).filter(
                                Interactor.uniprotkb == row['#ID(s) interactor A'].split(':')[1]).one())
                    if (session.query(Interactor).filter(
                            Interactor.id == row['ID(s) interactor B'].split(':')[1]).first() != None):
                        interactors.append(session.query(Interactor).filter(
                            Interactor.id == row['ID(s) interactor B'].split(':')[1]).one())
                    elif (session.query(Interactor).filter(
                                Interactor.uniprotkb == row['ID(s) interactor B'].split(':')[1]).first() != None):
                            interactors.append(session.query(Interactor).filter(
                                Interactor.uniprotkb == row['ID(s) interactor B'].split(':')[1]).one())

                    if (len(interactors) == 2):
                        interaction = session.query(Interaction).filter(
                            (Interaction.interactors.contains(interactors[0])),
                            (Interaction.interactors.contains(
                                interactors[1]))).first()
                        if (interaction == None):
                            interaction = Interaction(strain=strain,
                                                      type=(interactors[0].type + '-' + interactors[1].type))
                            if (interactors[0] == interactors[1]):
                                del interactors[1]
                            interaction.interactors = interactors
                            session.add(interaction)
                            session.commit()
                        reference = InteractionReference(interaction_type=row['Interaction type(s)'].split('(')[1][:-1],
                                                         interaction_id=interaction.id,
                                                         confidence_score=row['Confidence value(s)'],
                                                         pmid=row['Publication Identifier(s)'].split('med:')[1][:8],
                                                         detection_method=
                                                         row['Interaction detection method(s)'].split('(')[
                                                             1][
                                                         :-1],
                                                         author_last_name=row['Publication 1st author(s)'].split(',')[
                                                             0],
                                                         publication_date=row['Publication 1st author(s)'].split('(')[
                                                                              1][
                                                                          :-1],
                                                         publication_ref=row['Publication 1st author(s)'],
                                                         source_db=row['Source database(s)'].split('(')[1][:-1],
                                                         experimental_role_a=
                                                         row['Experimental role(s) interactor A'].split('(')[1][:-1],
                                                         experimental_role_b=
                                                         row['Experimental role(s) interactor B'].split('(')[1][:-1],
                                                         interactor_a_id=row['#ID(s) interactor A'].split(':')[1],
                                                         interactor_b_id=row['ID(s) interactor B'].split(':')[1])
                        session.add(reference)
                        interaction.is_experimental = 1
            session.commit()
            print(session.query(Interaction).count())


    def parse_IntAct(file, strain, taxid):
        with open(file) as csvfile:
            reader = csv.DictReader(csvfile, delimiter='\t')
            for row in reader:
                interactors = []

                if ((row['Taxid interactor A'].split('|')[0] == taxid) &
                        (row['Taxid interactor B'].split('|')[0] == taxid)):

                    if (session.query(Interactor).filter(
                            Interactor.id == row['#ID(s) interactor A'].split(':')[1]).first() != None):
                        interactors.append(session.query(Interactor).filter(
                            Interactor.id == row['#ID(s) interactor A'].split(':')[1]).one())
                    elif (session.query(Interactor).filter(
                                Interactor.uniprotkb == row['#ID(s) interactor A'].split(':')[1]).first() != None):
                            interactors.append(session.query(Interactor).filter(
                                Interactor.uniprotkb == row['#ID(s) interactor A'].split(':')[1]).one())
                    if (session.query(Interactor).filter(
                            Interactor.id == row['ID(s) interactor B'].split(':')[1]).first() != None):
                        interactors.append(session.query(Interactor).filter(
                            Interactor.id == row['ID(s) interactor B'].split(':')[1]).one())
                    elif (session.query(Interactor).filter(
                                Interactor.uniprotkb == row['ID(s) interactor B'].split(':')[1]).first() != None):
                            interactors.append(session.query(Interactor).filter(
                                Interactor.uniprotkb == row['ID(s) interactor B'].split(':')[1]).one())

                    if (len(interactors) == 2):
                        interaction = session.query(Interaction).filter(
                            (Interaction.interactors.contains(interactors[0])),
                            (Interaction.interactors.contains(
                                interactors[1]))).first()
                        if (interaction == None):
                            interaction = Interaction(strain=strain,
                                                      type=(interactors[0].type + '-' + interactors[1].type))
                            if (interactors[0] == interactors[1]):
                                del interactors[1]
                            interaction.interactors = interactors
                            session.add(interaction)
                            session.commit()
                        reference = InteractionReference(interaction_type=row['Interaction type(s)'].split('(')[1][:-1],
                                                         interaction_id=interaction.id,
                                                         confidence_score=row['Confidence value(s)'],
                                                         pmid=row['Publication Identifier(s)'].split('med:')[1][:8],
                                                         detection_method=
                                                         row['Interaction detection method(s)'].split('(')[
                                                             1][
                                                         :-1],
                                                         author_last_name=row['Publication 1st author(s)'].split(',')[
                                                             0],
                                                         publication_date=row['Publication 1st author(s)'].split('(')[
                                                                              1][
                                                                          :-1],
                                                         publication_ref=row['Publication 1st author(s)'],
                                                         source_db=row['Source database(s)'].split('(')[1][:-1],
                                                         experimental_role_a=
                                                         row['Experimental role(s) interactor A'].split('(')[1][:-1],
                                                         experimental_role_b=
                                                         row['Experimental role(s) interactor B'].split('(')[1][:-1],
                                                         interactor_a_id=row['#ID(s) interactor A'].split(':')[1],
                                                         interactor_b_id=row['ID(s) interactor B'].split(':')[1])
                        session.add(reference)
                        interaction.is_experimental = 1
            session.commit()
            print(session.query(Interaction).count())


    def parse_iRefIndex(file, strain, taxid):
        with open(file) as csvfile:
            reader = csv.DictReader(csvfile, delimiter='\t')
            for row in reader:
                interactors = []

                if ((row['Taxid interactor A'].split('|')[0] == taxid) &
                        (row['Taxid interactor B'].split('|')[0] == taxid)):

                    if (row['#ID(s) interactor A'].split(':')[0] == 'uniprotkb'):
                        if (session.query(Interactor).filter(
                                Interactor.id == row['#ID(s) interactor A'].split(':')[1]).first() != None):
                            interactors.append(session.query(Interactor).filter(
                                Interactor.id == row['#ID(s) interactor A'].split(':')[1]).one())
                        else:
                            if (session.query(Interactor).filter(
                                    Interactor.uniprotkb == row['#ID(s) interactor A'].split(':')[1]).first() != None):
                                interactors.append(session.query(Interactor).filter(
                                    Interactor.uniprotkb == row['#ID(s) interactor A'].split(':')[1]).one())
                    if (row['#ID(s) interactor A'].split(':')[0] == 'refseq'):
                        if (session.query(InteractorXref).filter(InteractorXref.accession ==
                                                                 row['#ID(s) interactor A'].split(':')[
                                                                     1]).first() != None):
                            interactors.append(session.query(InteractorXref).filter(InteractorXref.accession ==
                                                                                    row['#ID(s) interactor A'].split(
                                                                                        ':')[
                                                                                        1]).one().interactor)
                    if (row['ID(s) interactor B'].split(':')[0] == 'uniprotkb'):
                        if (session.query(Interactor).filter(
                                Interactor.id == row['ID(s) interactor B'].split(':')[1]).first() != None):
                            interactors.append(session.query(Interactor).filter(
                                Interactor.id == row['ID(s) interactor B'].split(':')[1]).one())
                        else:
                            if (session.query(Interactor).filter(
                                    Interactor.uniprotkb == row['ID(s) interactor B'].split(':')[1]).first() != None):
                                interactors.append(session.query(Interactor).filter(
                                    Interactor.uniprotkb == row['ID(s) interactor B'].split(':')[1]).one())
                    if (row['ID(s) interactor B'].split(':')[0] == 'refseq'):
                        if (session.query(InteractorXref).filter(InteractorXref.accession ==
                                                                 row['ID(s) interactor B'].split(':')[
                                                                     1]).first() != None):
                            interactors.append(session.query(InteractorXref).filter(InteractorXref.accession ==
                                                                                    row['ID(s) interactor B'].split(
                                                                                        ':')[1]).one().interactor)

                    if (len(interactors) == 2):
                        interaction = session.query(Interaction).filter(
                            (Interaction.interactors.contains(interactors[0])),
                            (Interaction.interactors.contains(
                                interactors[1]))).first()
                        if (interaction == None):
                            interaction = Interaction(strain=strain,
                                                      type=(interactors[0].type + '-' + interactors[1].type))
                            if (interactors[0] == interactors[1]):
                                del interactors[1]
                            interaction.interactors = interactors
                            session.add(interaction)
                            session.commit()
                        author, date, ref, type, pmid, detection = None, None, None, None, None, None
                        if (row['Publication 1st author(s)'] != '-'):
                            author = row['Publication 1st author(s)'].split('-')[0]
                            date = row['Publication 1st author(s)'].split('-')[1]
                            ref = row['Publication 1st author(s)']
                        if (row['Interaction type(s)'] != '-'):
                            type = row['Interaction type(s)'].split('(')[1][:-1]
                        if (row['Publication Identifier(s)'] != '-'):
                            pmid = row['Publication Identifier(s)'].split('med:')[1][:8]
                        if (row['Interaction detection method(s)'] != '-'):
                            detection = row['Interaction detection method(s)'].split('(')[1][:-1]
                        # there are more than one pmid sometimes
                        reference = InteractionReference(interaction_type=type,
                                                         interaction_id=interaction.id,
                                                         confidence_score=row['Confidence value(s)'],
                                                         pmid=pmid,
                                                         detection_method=detection,
                                                         author_last_name=author,
                                                         publication_date=date,
                                                         publication_ref=ref,
                                                         source_db=row['Source database(s)'].split('(')[1][:-1],
                                                         experimental_role_a=
                                                         row['Experimental role(s) interactor A'].split('(')[1][:-1],
                                                         experimental_role_b=
                                                         row['Experimental role(s) interactor B'].split('(')[1][:-1],
                                                         interactor_a_id=row['#ID(s) interactor A'].split(':')[1],
                                                         interactor_b_id=row['ID(s) interactor B'].split(':')[1])
                        session.add(reference)
                        interaction.is_experimental = 1
            session.commit()
            print(session.query(Interaction).count())


    def parse_mentha(file, strain, taxid):
        with open(file) as csvfile:
            reader = csv.DictReader(csvfile, delimiter='\t')
            for row in reader:
                interactors = []

                if ((row['Taxid interactor A'].split('|')[0] == taxid) &
                        (row['Taxid interactor B'].split('|')[0] == taxid)):

                    if (session.query(Interactor).filter(
                            Interactor.id == row['#ID(s) interactor A'].split(':')[1]).first() != None):
                        interactors.append(session.query(Interactor).filter(
                            Interactor.id == row['#ID(s) interactor A'].split(':')[1]).one())
                    else:
                        if (session.query(Interactor).filter(
                                Interactor.uniprotkb == row['#ID(s) interactor A'].split(':')[1]).first() != None):
                            interactors.append(session.query(Interactor).filter(
                                Interactor.uniprotkb == row['#ID(s) interactor A'].split(':')[1]).one())
                    if (session.query(Interactor).filter(
                            Interactor.id == row['ID(s) interactor B'].split(':')[1]).first() != None):
                        interactors.append(session.query(Interactor).filter(
                            Interactor.id == row['ID(s) interactor B'].split(':')[1]).one())
                    else:
                        if (session.query(Interactor).filter(
                                Interactor.uniprotkb == row['ID(s) interactor B'].split(':')[1]).first() != None):
                            interactors.append(session.query(Interactor).filter(
                                Interactor.uniprotkb == row['ID(s) interactor B'].split(':')[1]).one())

                    if (len(interactors) == 2):
                        interaction = session.query(Interaction).filter(
                            (Interaction.interactors.contains(interactors[0])),
                            (Interaction.interactors.contains(
                                interactors[1]))).first()
                        if (interaction == None):
                            interaction = Interaction(strain=strain,
                                                      type=(interactors[0].type + '-' + interactors[1].type))
                            if (interactors[0] == interactors[1]):
                                del interactors[1]
                            interaction.interactors = interactors
                            session.add(interaction)
                            session.commit()
                        reference = InteractionReference(interaction_type=row['Interaction type(s)'].split('(')[1][:-1],
                                                         interaction_id=interaction.id,
                                                         confidence_score=row['Confidence value(s)'],
                                                         pmid=row['Publication Identifier(s)'].split('med:')[1][:8],
                                                         detection_method=
                                                         row['Interaction detection method(s)'].split('(')[1][:-1],
                                                         source_db=row['Source database(s)'].split('(')[1][:-1])
                        session.add(reference)
                        interaction.is_experimental = 1
            session.commit()
            print(session.query(Interaction).count())


    def parse_MINT(file, strain, taxid):
        with open(file) as csvfile:
            reader = csv.DictReader(csvfile, delimiter='\t')
            for row in reader:
                interactors = []

                if ((row['Taxid interactor A'].split('|')[0] == taxid) &
                        (row['Taxid interactor B'].split('|')[0] == taxid)):

                    if (session.query(Interactor).filter(
                            Interactor.id == row['#ID(s) interactor A'].split(':')[1]).first() != None):
                        interactors.append(session.query(Interactor).filter(
                            Interactor.id == row['#ID(s) interactor A'].split(':')[1]).one())
                    else:
                        if (session.query(Interactor).filter(
                                Interactor.uniprotkb == row['#ID(s) interactor A'].split(':')[1]).first() != None):
                            interactors.append(session.query(Interactor).filter(
                                Interactor.uniprotkb == row['#ID(s) interactor A'].split(':')[1]).one())
                    if (session.query(Interactor).filter(
                            Interactor.id == row['ID(s) interactor B'].split(':')[1]).first() != None):
                        interactors.append(session.query(Interactor).filter(
                            Interactor.id == row['ID(s) interactor B'].split(':')[1]).one())
                    else:
                        if (session.query(Interactor).filter(
                                Interactor.uniprotkb == row['ID(s) interactor B'].split(':')[1]).first() != None):
                            interactors.append(session.query(Interactor).filter(
                                Interactor.uniprotkb == row['ID(s) interactor B'].split(':')[1]).one())

                    if (len(interactors) == 2):
                        interaction = session.query(Interaction).filter(
                            (Interaction.interactors.contains(interactors[0])),
                            (Interaction.interactors.contains(
                                interactors[1]))).first()
                        if (interaction == None):
                            interaction = Interaction(strain=strain,
                                                      type=(interactors[0].type + '-' + interactors[1].type))
                            if (interactors[0] == interactors[1]):
                                del interactors[1]
                            interaction.interactors = interactors
                            session.add(interaction)
                            session.commit()
                        reference = InteractionReference(interaction_type=row['Interaction type(s)'].split('(')[1][:-1],
                                                         interaction_id=interaction.id,
                                                         confidence_score=row['Confidence value(s)'],
                                                         pmid=row['Publication Identifier(s)'].split('med:')[1][:8],
                                                         detection_method=
                                                         row['Interaction detection method(s)'].split('(')[
                                                             1][
                                                         :-1],
                                                         author_last_name=row['Publication 1st author(s)'].split(',')[
                                                             0],
                                                         publication_date=row['Publication 1st author(s)'].split('(')[
                                                                              1][
                                                                          :-1],
                                                         publication_ref=row['Publication 1st author(s)'],
                                                         source_db=row['Source database(s)'].split('(')[1][:-1],
                                                         experimental_role_a=
                                                         row['Experimental role(s) interactor A'].split('(')[1][:-1],
                                                         experimental_role_b=
                                                         row['Experimental role(s) interactor B'].split('(')[1][:-1],
                                                         interactor_a_id=row['#ID(s) interactor A'].split(':')[1],
                                                         interactor_b_id=row['ID(s) interactor B'].split(':')[1])
                        session.add(reference)
                        interaction.is_experimental = 1
            session.commit()
            print(session.query(Interaction).count())

    parse_IMEx('PSICQUIC_Paeruginosa_PAO1_IMEx.txt', 'PAO1', 'taxid:208964(pseae)')
    parse_IMEx('PSICQUIC_Paeruginosa_PA14_IMEx.txt', 'PA14', 'taxid:208963(pseab)')

    parse_iRefIndex('PSICQUIC_Paeruginosa_PAO1_iRefIndex.txt', 'PAO1', 'taxid:208964(Pseudomonas aeruginosa PAO1)')
    parse_iRefIndex('PSICQUIC_Paeruginosa_PA14_iRefIndex.txt', 'PA14',
                    'taxid:208963(Pseudomonas aeruginosa UCBPP-PA14)')

    parse_mentha('PSICQUIC_Paeruginosa_PAO1_mentha.txt', 'PAO1', 'taxid:208964(Pseudomonas aeruginosa PAO1)')
    parse_mentha('PSICQUIC_Paeruginosa_PA14_mentha.txt', 'PA14', 'taxid:208963(Pseudomonas aeruginosa UCBPP-PA14)')

    parse_MINT('PSICQUIC_Paeruginosa_PAO1_MINT.txt', 'PAO1', 'taxid:208964(pseae)')
    parse_MINT('PSICQUIC_Paeruginosa_PA14_MINT.txt', 'PA14', 'taxid:208963(pseab)')

    parse_IntAct('PSICQUIC_Paeruginosa_PAO1_IntAct.txt', 'PAO1', 'taxid:208964(pseae)')
    parse_IntAct('IntAct_Paeruginosa_PA14.txt', 'PA14', 'taxid:208963(pseab)')
    parse_IntAct('IntAct_Paeruginosa_PAO1.txt', 'PAO1', 'taxid:208964(pseae)')

    with open('interactions_in_P_aeruginosa_geoff.csv') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            interactors = []
            if (session.query(Interactor).filter(Interactor.id == row['locus_tag']).first() != None):
                interactors.append(session.query(Interactor).filter(Interactor.id == row['locus_tag']).one())
            locus_tag1 = row['locus_tag']
            row = next(reader)
            if (session.query(Interactor).filter(Interactor.id == row['locus_tag']).first() != None):
                interactors.append(session.query(Interactor).filter(Interactor.id == row['locus_tag']).one())

            if (len(interactors) == 2):
                interaction = session.query(Interaction).filter((Interaction.interactors.contains(interactors[0])),
                                                                (Interaction.interactors.contains(
                                                                    interactors[1]))).first()
                if (interaction == None):
                    interaction = Interaction(strain='PAO1', type=(interactors[0].type + '-' + interactors[1].type))
                    if (interactors[0] == interactors[1]):
                        del interactors[1]
                    interaction.interactors = interactors
                    session.add(interaction)
                    session.commit()

                reference = InteractionReference(interaction_type=row['type'], interaction_id=interaction.id,
                                                 pmid=row['pmid'], interaction_full_name=row['full_name'],
                                                 detection_method=row['experimental_type'])
                interaction.is_experimental = 1
                session.add(reference)
        session.commit()
        print(session.query(Interaction).count())

    with open('208964.psicquic-mitab_2.5.v10.5.txt') as csvfile:
        fieldnames = ['interactor_A', 'interactor_B', 'altID_A', 'altID_B', 'alias_A', 'alias_B', 'detection',
                      'publication', 'publication_ID', 'taxid_A', 'taxid_B', 'type', 'source_db', 'identifier',
                      'confidence']
        reader = csv.DictReader(csvfile, delimiter='\t', fieldnames=fieldnames)
        for row in reader:
            interactors = []
            locus_tag1 = row['interactor_A'].split('|')[0].split('.')[1]
            if (session.query(Interactor).filter(Interactor.id == locus_tag1).first() != None):
                interactors.append(session.query(Interactor).filter(Interactor.id == locus_tag1).one())
            locus_tag2 = row['interactor_B'].split('|')[0].split('.')[1]
            if (session.query(Interactor).filter(Interactor.id == locus_tag2).first() != None):
                interactors.append(session.query(Interactor).filter(Interactor.id == locus_tag2).one())

            if (len(interactors) == 2):
                interaction = session.query(Interaction).filter((Interaction.interactors.contains(interactors[0])),
                                                                (Interaction.interactors.contains(
                                                                    interactors[1]))).first()
                if (interaction == None):
                    interaction = Interaction(strain='PAO1', type=(interactors[0].type + '-' + interactors[1].type),
                                              is_experimental=0)
                    if (interactors[0] == interactors[1]):
                        del interactors[1]
                    interaction.interactors = interactors
                    session.add(interaction)
                    session.commit()

                detection = row['detection'].split('(')[1][:-1]
                confidence = row['confidence'].split(':')[1]
                type = find_type(row['type'].split(':')[2][:-1])
                source_db, author, date, ref, pmid = None, None, None, None, None

                if (row['source_db'] != '-'):
                    source_db = row['source_db'].split('(')[1][:-1]
                if (row['publication'] != '-'):
                    author_ln = row['publication'].split(' ')[0]
                    publication_date = row['publication'].split('(')[1][:-1]
                    publication_ref = row['publication']
                if (row['publication_ID'] != '-'):
                    pmid = row['publication_ID'].split(':')[1]
                reference = InteractionReference(interaction_type=type, interaction_id=interaction.id,
                                                 detection_method=detection, confidence_score=confidence,
                                                 source_db=source_db, pmid=pmid, author_last_name=author,
                                                 publication_ref=ref, publication_date=date)

                session.add(reference)
                if (detection == 'experimental interaction detection'):
                    interaction.is_experimental = 1
        session.commit()
        print(session.query(Interaction).count())

    with open('PAO1_Interactome_predicted.csv') as csvfile:
        reader = csv.DictReader(csvfile)

        for row in reader:
            if (float(row['Confidence']) >= 0.9):
                interactors = []
                if (session.query(Interactor).filter(Interactor.id == row['Protein1']).first() != None):
                    interactors.append(session.query(Interactor).filter(Interactor.id == row['Protein1']).one())

                if (session.query(Interactor).filter(Interactor.id == row['Protein2']).first() != None):
                    interactors.append(session.query(Interactor).filter(Interactor.id == row['Protein2']).one())

                if (len(interactors) == 2):
                    interaction = session.query(Interaction).filter((Interaction.interactors.contains(interactors[0])),
                                                                    (Interaction.interactors.contains(
                                                                        interactors[1]))).first()
                    if (interaction == None):
                        interaction = Interaction(strain='PAO1', type=(interactors[0].type + '-' + interactors[1].type),
                                                  is_experimental=0)
                        if (interactors[0] == interactors[1]):
                            del interactors[1]
                        interaction.interactors = interactors
                        session.add(interaction)
                        session.commit()

                    reference = InteractionReference(interaction_type='predicted', interaction_id=interaction.id,
                                                     confidence_score=row['Confidence'], pmid='22848443',
                                                     detection_method='computational prediction')
                    session.add(reference)
        session.commit()
        print(session.query(Interaction).count())

    with open('xlinkdb-ouput_PAO1.txt') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
            interactors = []
            if (session.query(Interactor).filter(Interactor.id == row['proA']).first() != None):
                interactors.append(session.query(Interactor).filter(Interactor.id == row['proA']).one())
            else:
                if (session.query(Interactor).filter(Interactor.uniprotkb == row['proA']).first() != None):
                    interactors.append(session.query(Interactor).filter(Interactor.uniprotkb == row['proA']).one())
            if (session.query(Interactor).filter(Interactor.id == row['proB']).first() != None):
                interactors.append(session.query(Interactor).filter(Interactor.id == row['proB']).one())
            else:
                if (session.query(Interactor).filter(Interactor.uniprotkb == row['proB']).first() != None):
                    interactors.append(session.query(Interactor).filter(Interactor.uniprotkb == row['proB']).one())

            if (len(interactors) == 2):
                interaction = session.query(Interaction).filter((Interaction.interactors.contains(interactors[0])),
                                                                (Interaction.interactors.contains(
                                                                    interactors[1]))).first()
                if (interaction == None):
                    interaction = Interaction(strain='PAO1', type=(interactors[0].type + '-' + interactors[1].type))
                    if (interactors[0] == interactors[1]):
                        del interactors[1]
                    interaction.interactors = interactors
                    session.add(interaction)
                    session.commit()
                reference = InteractionReference(interaction_type='physical association', interaction_id=interaction.id,
                                                 pmid='25800553', source_db='XLinkDB',
                                                 detection_method='chemical cross-linking mass spectrometry')
                session.add(reference)
                interaction.is_experimental = 1
        session.commit()
        print(session.query(Interaction).count())

    with open('PSICQUIC_Paeruginosa_PAO1_ADIPInteractomes.txt') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
            interactors = []

            if ((row['Taxid interactor A'] == 'taxid:208964(Pseudomonas aeruginosa PAO1)') &
                    (row['Taxid interactor B'] == 'taxid:208964(Pseudomonas aeruginosa PAO1)')):

                if (session.query(Interactor).filter(
                        Interactor.id == row['#ID(s) interactor A'].split(':')[1]).first() != None):
                    interactors.append(session.query(Interactor).filter(
                        Interactor.id == row['#ID(s) interactor A'].split(':')[1]).one())
                else:
                    if (session.query(Interactor).filter(
                            Interactor.uniprotkb == row['#ID(s) interactor A'].split(':')[1]).first() != None):
                        interactors.append(session.query(Interactor).filter(
                            Interactor.uniprotkb == row['#ID(s) interactor A'].split(':')[1]).one())

                if (session.query(Interactor).filter(
                        Interactor.id == row['ID(s) interactor B'].split(':')[1]).first() != None):
                    interactors.append(session.query(Interactor).filter(
                        Interactor.id == row['ID(s) interactor B'].split(':')[1]).one())
                else:
                    if (session.query(Interactor).filter(
                            Interactor.uniprotkb == row['ID(s) interactor B'].split(':')[1]).first() != None):
                        interactors.append(session.query(Interactor).filter(
                            Interactor.uniprotkb == row['ID(s) interactor B'].split(':')[1]).one())

                if (len(interactors) == 2):
                    interaction = session.query(Interaction).filter((Interaction.interactors.contains(interactors[0])),
                                                                    (Interaction.interactors.contains(
                                                                        interactors[1]))).first()
                    if (interaction == None):
                        interaction = Interaction(strain='PAO1',
                                                  type=(interactors[0].typ33e + '-' + interactors[1].type))
                        if (interactors[0] == interactors[1]):
                            del interactors[1]
                        interaction.interactors = interactors
                        session.add(interaction)
                        session.commit()
                    reference = InteractionReference(interaction_type=row['Interaction type(s)'],
                                                     interaction_id=interaction.id,
                                                     confidence_score=row['Confidence value(s)'],
                                                     pmid=row['Publication Identifier(s)'].split('med:')[1][:8],
                                                     detection_method=row['Interaction detection method(s)'].split('(')[
                                                                          1][:-1],
                                                     author_last_name=row['Publication 1st author(s)'].split(',')[0],
                                                     publication_date=row['Publication 1st author(s)'].split('(')[1][
                                                                      :-1],
                                                     publication_ref=row['Publication 1st author(s)'])
                    session.add(reference)
                    interaction.is_experimental = 1
                    xref = session.query(InteractionXref).filter(
                        InteractionXref.accession == row['Interaction identifier(s)'].split(':')[1]).first()
                    if (xref == None):
                        xref = InteractionXref(accession=row['Interaction identifier(s)'].split(':')[1],
                                               data_source=row['Source database(s)'].split('(')[1][:-1],
                                               interaction_id=interaction.id)
                        session.add(xref)
        session.commit()
        print(session.query(Interaction).count())

    with open('PSICQUIC_Paeruginosa_PAO1_MPIDB.txt') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
            interactors = []

            if ((row['Taxid interactor A'].split('|')[0] == 'taxid:208964(pseae)') &
                    (row['Taxid interactor B'].split('|')[0] == 'taxid:208964(pseae)')):

                if (session.query(Interactor).filter(
                        Interactor.id == row['#ID(s) interactor A'].split(':')[1]).first() != None):
                    interactors.append(session.query(Interactor).filter(
                        Interactor.id == row['#ID(s) interactor A'].split(':')[1]).one())
                else:
                    if (session.query(Interactor).filter(
                            Interactor.uniprotkb == row['#ID(s) interactor A'].split(':')[1]).first() != None):
                        interactors.append(session.query(Interactor).filter(
                            Interactor.uniprotkb == row['#ID(s) interactor A'].split(':')[1]).one())
                if (session.query(Interactor).filter(
                        Interactor.id == row['ID(s) interactor B'].split(':')[1]).first() != None):
                    interactors.append(session.query(Interactor).filter(
                        Interactor.id == row['ID(s) interactor B'].split(':')[1]).one())
                else:
                    if (session.query(Interactor).filter(
                            Interactor.uniprotkb == row['ID(s) interactor B'].split(':')[1]).first() != None):
                        interactors.append(session.query(Interactor).filter(
                            Interactor.uniprotkb == row['ID(s) interactor B'].split(':')[1]).one())

                if (len(interactors) == 2):
                    interaction = session.query(Interaction).filter((Interaction.interactors.contains(interactors[0])),
                                                                    (Interaction.interactors.contains(
                                                                        interactors[1]))).first()
                    if (interaction == None):
                        interaction = Interaction(strain='PAO1', type=(interactors[0].type + '-' + interactors[1].type))
                        if (interactors[0] == interactors[1]):
                            del interactors[1]
                        interaction.interactors = interactors
                        session.add(interaction)
                        session.commit()
                    reference = InteractionReference(interaction_type=row['Interaction type(s)'].split('(')[1][:-1],
                                                     interaction_id=interaction.id,
                                                     confidence_score=row['Confidence value(s)'],
                                                     pmid=row['Publication Identifier(s)'].split('med:')[1][:8],
                                                     detection_method=row['Interaction detection method(s)'].split('(')[
                                                                          1][
                                                                      :-1],
                                                     author_last_name=row['Publication 1st author(s)'].split(' ')[0],
                                                     publication_date=row['Publication 1st author(s)'].split('(')[1][
                                                                      :-1],
                                                     publication_ref=row['Publication 1st author(s)'],
                                                     source_db=row['Source database(s)'].split('(')[1][:-1],
                                                     experimental_role_a=
                                                     row['Experimental role(s) interactor A'].split('(')[1][:-1],
                                                     experimental_role_b=
                                                     row['Experimental role(s) interactor B'].split('(')[1][:-1],
                                                     interactor_a_id=row['#ID(s) interactor A'].split(':')[1],
                                                     interactor_b_id=row['ID(s) interactor B'].split(':')[1])
                    session.add(reference)
                    interaction.is_experimental = 1
        session.commit()
        print(session.query(Interaction).count())

    with open('Paeruginosa_regulatory_network.csv') as csvfile:
        reader = csv.DictReader(csvfile)
        string = []
        for row in reader:
            strains = row['Strain'].split(',')
            for strain in strains:
                if ((strain == 'PA14') | (strain == 'PAO1')):
                    interactors = []
                    if (session.query(Interactor).filter((Interactor.name == row['Regulator (TF or sigma)']),
                                                         (Interactor.strain == strain)).count() == 1):
                        interactors.append(
                            session.query(Interactor).filter((Interactor.name == row['Regulator (TF or sigma)']),
                                                             (Interactor.strain == strain)).one())
                    else:
                        if (session.query(Interactor).filter(
                                Interactor.id == row['Regulator (TF or sigma)']).count() == 1):
                            interactors.append(
                                session.query(Interactor).filter(Interactor.id == row['Regulator (TF or sigma)']).one())
                    if (session.query(Interactor).filter((Interactor.name == row['Target']),
                                                         (Interactor.strain == strain)).count() == 1):
                        interactors.append(session.query(Interactor).filter((Interactor.name == row['Target']),
                                                                            (Interactor.strain == strain)).one())
                    else:
                        if (session.query(Interactor).filter(Interactor.id == row['Target']).count() == 1):
                            interactors.append(
                                session.query(Interactor).filter(Interactor.id == row['Target']).one())

                    string.append(len(interactors))
                    if (len(interactors) == 2):
                        # string.append(row['Regulator (TF or sigma)'] + ':' + row['Target'] + '(' + strain + ')')
                        interaction = session.query(Interaction).filter(
                            (Interaction.interactors.contains(interactors[0])),
                            (Interaction.interactors.contains(
                                interactors[1]))).first()
                        if (interactors[0].type != 'CDS(TF/sigma)'):
                            interactors[0].type = 'CDS(TF/sigma)'

                        if (interaction != None):
                            if (interaction.type != 'TF/sigma-binding site'):
                                interaction.type = 'TF/sigma-binding site'
                                session.commit()

                        second_interactor = interactors[1]
                        if (interaction == None):
                            interaction = Interaction(strain=strain, type='TF/sigma-binding site')
                            if (interactors[0] == interactors[1]):
                                del interactors[1]
                            interaction.interactors = interactors
                            session.add(interaction)
                            session.commit()

                        reference = InteractionReference(interaction_id=interaction.id, pmid=row['pmid'],
                                                         interaction_type='TF/sigma-binding site (' + row['mode'] +
                                                                          'regulation)',
                                                         detection_method=row['evidence'],
                                                         interaction_full_name=interactors[0].id + ' regulates(' +
                                                                               row[
                                                                                   'mode'] + ') ' + second_interactor.id,
                                                         source_db=row['source_db'])
                        session.add(reference)
                        interaction.is_experimental = 1
        session.commit()
        print(session.query(InteractionReference).count())

    with open('PAO1-PA14_orthologs.csv') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if (session.query(Interactor).filter(Interactor.id == row['Locus Tag (Strain 1)']).first() != None):
                if (session.query(Interactor).filter(Interactor.id == row['Locus Tag (Strain 2)']).first() != None):
                    ortholog1 = OrthologPseudomonas(interactor_id=row['Locus Tag (Strain 1)'], strain_interactor='PAO1',
                                         ortholog_id=row['Locus Tag (Strain 2)'], strain_ortholog='PA14')
                    ortholog2 = OrthologPseudomonas(interactor_id=row['Locus Tag (Strain 2)'], strain_interactor='PA14',
                                         ortholog_id=row['Locus Tag (Strain 1)'], strain_ortholog='PAO1')

                    session.add(ortholog1)
                    session.add(ortholog2)
        session.commit()

    def orthologs_ecoli(file, strain):
        with open(file) as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                if (session.query(Interactor).filter(Interactor.id == row['Locus Tag (Strain 2)']).first() != None):
                    ortholog = OrthologEcoli(interactor_id=row['Locus Tag (Strain 2)'], strain_interactor=strain,
                                        ortholog_id=row['Locus Tag (Strain 1)'], strain_ortholog='E. coli K12',
                                             ortholog_refseq = row['NCBI RefSeq Accession (Strain 1)'])
                    if (row['Ortholuge Class'] == ''):
                        ortholog.ortholuge_classification = 'RBBH'
                    else:
                        ortholog.ortholuge_classification = row['Ortholuge Class']
                    session.add(ortholog)
            session.commit()

    orthologs_ecoli('PAO1-EcoliK12MG1655_orthologs.csv', 'PAO1')
    orthologs_ecoli('PA14-EcoliK12MG1655_orthologs.csv', 'PA14')

    with open('Ecoli_IDs.csv') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if (row['Gene names'] != ''):
                locus = row['Gene names'].split(' ')[0]
                strains = ['PAO1', 'PA14']
                for strain in strains:
                    if (session.query(OrthologEcoli).filter((OrthologEcoli.ortholog_id == locus),
                                                            (OrthologEcoli.strain_interactor == strain)).first() != None):
                        ortholog = session.query(OrthologEcoli).filter((OrthologEcoli.ortholog_id == locus),
                                                            (OrthologEcoli.strain_interactor == strain)).one()
                        ortholog.ortholog_uniprot = row['Entry']
                        ortholog.name = row['Gene name']
                session.commit()

    def parse_Ecoli_IntAct(file):
        with open(file) as csvfile:
            reader = csv.DictReader(csvfile, delimiter='\t')
            for row in reader:
                interactorsPAO1, interactorsPA14 = [], []
                interactors = [interactorsPAO1, interactorsPA14]
                strains = ['PAO1', 'PA14']

                for num in range(0,2):
                    if (session.query(OrthologEcoli).filter(
                            (OrthologEcoli.ortholog_uniprot == row['#ID(s) interactor A'].split(':')[1]),
                            (OrthologEcoli.strain_interactor == strains[num])).first() != None):
                        id = session.query(OrthologEcoli).filter(
                            (OrthologEcoli.ortholog_uniprot == row['#ID(s) interactor A'].split(':')[1]),
                            (OrthologEcoli.strain_interactor == strains[num])).one().interactor_id
                        interactors[num].append(session.query(Interactor).filter(Interactor.id == id).one())
                    if (session.query(OrthologEcoli).filter(
                            (OrthologEcoli.ortholog_uniprot == row['ID(s) interactor B'].split(':')[1]),
                            (OrthologEcoli.strain_interactor == strains[num])).first() != None):
                        id = session.query(OrthologEcoli).filter(
                            (OrthologEcoli.ortholog_uniprot == row['ID(s) interactor B'].split(':')[1]),
                            (OrthologEcoli.strain_interactor == strains[num])).one().interactor_id
                        interactors[num].append(session.query(Interactor).filter(Interactor.id == id).one())

                    if (len(interactors[num]) == 2):
                        interaction = session.query(Interaction).filter((Interaction.interactors.contains(interactors[num][0])),
                                                                        (Interaction.interactors.contains(
                                                                            interactors[num][1]))).first()
                        if (interaction != None):
                            interaction.is_Ecoli_ortholog = 2
                        if (interaction == None):
                            interaction = Interaction(strain=strains[num], type=(interactors[num][0].type + '-' + interactors[num][
                                1].type), is_experimental=0, is_Ecoli_ortholog=1)
                            if (interactors[num][0] == interactors[num][1]):
                                del interactors[num][1]
                            interaction.interactors = interactors[num]
                            session.add(interaction)
                            session.commit()
                        reference = InteractionReference(interaction_type=row['Interaction type(s)'].split('(')[1][:-1],
                                                         interaction_id=interaction.id,
                                                         confidence_score=row['Confidence value(s)'],
                                                         pmid=row['Publication Identifier(s)'].split('med:')[1][:8],
                                                         detection_method=row['Interaction detection method(s)'].split('(')[
                                                                              1][
                                                                          :-1],
                                                         author_last_name=row['Publication 1st author(s)'].split(',')[0],
                                                         publication_date=row['Publication 1st author(s)'].split('(')[1][
                                                                          :-1],
                                                         publication_ref=row['Publication 1st author(s)'],
                                                         source_db=row['Source database(s)'].split('(')[1][:-1],
                                                         experimental_role_a=
                                                         row['Experimental role(s) interactor A'].split('(')[1][:-1],
                                                         experimental_role_b=
                                                         row['Experimental role(s) interactor B'].split('(')[1][:-1],
                                                         interactor_a_id=row['#ID(s) interactor A'].split(':')[1],
                                                         interactor_b_id=row['ID(s) interactor B'].split(':')[1])
                        session.add(reference)
            session.commit()
            print(session.query(Interaction).count())

    with open('PSICQUIC_Ecoli_IMEx.txt') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
            interactorsPAO1, interactorsPA14 = [], []
            interactors = [interactorsPAO1, interactorsPA14]
            strains = ['PAO1', 'PA14']

            for num in range(0, 2):
                if (session.query(OrthologEcoli).filter(
                        (OrthologEcoli.ortholog_uniprot == row['#ID(s) interactor A'].split(':')[1]),
                        (OrthologEcoli.strain_interactor == strains[num])).first() != None):
                    id = session.query(OrthologEcoli).filter(
                        (OrthologEcoli.ortholog_uniprot == row['#ID(s) interactor A'].split(':')[1]),
                        (OrthologEcoli.strain_interactor == strains[num])).one().interactor_id
                    interactors[num].append(session.query(Interactor).filter(Interactor.id == id).one())
                if (session.query(OrthologEcoli).filter(
                        (OrthologEcoli.ortholog_uniprot == row['ID(s) interactor B'].split(':')[1]),
                        (OrthologEcoli.strain_interactor == strains[num])).first() != None):
                    id = session.query(OrthologEcoli).filter(
                        (OrthologEcoli.ortholog_uniprot == row['ID(s) interactor B'].split(':')[1]),
                        (OrthologEcoli.strain_interactor == strains[num])).one().interactor_id
                    interactors[num].append(session.query(Interactor).filter(Interactor.id == id).one())

                if (len(interactors[num]) == 2):
                    interaction = session.query(Interaction).filter((Interaction.interactors.contains(interactors[num][0])),
                                                                    (Interaction.interactors.contains(
                                                                        interactors[num][1]))).first()
                    if (interaction != None):
                        interaction.is_Ecoli_ortholog = 2
                    if (interaction == None):
                        interaction = Interaction(strain=strains[num], type=(interactors[num][0].type + '-' + interactors[num][
                            1].type), is_experimental=0, is_Ecoli_ortholog = 1)
                        if (interactors[num][0] == interactors[num][1]):
                            del interactors[num][1]
                        interaction.interactors = interactors[num]
                        session.add(interaction)
                        session.commit()
                    reference = InteractionReference(interaction_type=row['Interaction type(s)'].split('(')[1][:-1],
                                                     interaction_id=interaction.id,
                                                     confidence_score=row['Confidence value(s)'],
                                                     pmid=row['Publication Identifier(s)'].split('med:')[1][:8],
                                                     detection_method=row['Interaction detection method(s)'].split('(')[
                                                                          1][:-1],
                                                     author_last_name=row['Publication 1st author(s)'].split(',')[0],
                                                     publication_date=row['Publication 1st author(s)'].split('(')[1][
                                                                      :-1],
                                                     publication_ref=row['Publication 1st author(s)'],
                                                     source_db=row['Source database(s)'].split('(')[1][:-1],
                                                     experimental_role_a=
                                                     row['Experimental role(s) interactor A'].split('(')[1][:-1],
                                                     experimental_role_b=
                                                     row['Experimental role(s) interactor B'].split('(')[1][:-1],
                                                     interactor_a_id=row['#ID(s) interactor A'].split(':')[1],
                                                     interactor_b_id=row['ID(s) interactor B'].split(':')[1])
                    session.add(reference)
        session.commit()
        print(session.query(Interaction).count())

    parse_Ecoli_IntAct('PSICQUIC_Ecoli_IntAct.txt')
    parse_Ecoli_IntAct('IntAct_Ecoli.txt')

    with open('PSICQUIC_Ecoli_iRefIndex.txt') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
            interactorsPAO1, interactorsPA14 = [], []
            interactors = [interactorsPAO1, interactorsPA14]
            strains = ['PAO1', 'PA14']

            for num in range(0, 2):
                if (row['#ID(s) interactor A'].split(':')[0] == 'uniprotkb'):
                    if (session.query(OrthologEcoli).filter(
                            (OrthologEcoli.ortholog_uniprot == row['#ID(s) interactor A'].split(':')[1]),
                            (OrthologEcoli.strain_interactor == strains[num])).first() != None):
                        id = session.query(OrthologEcoli).filter(
                            (OrthologEcoli.ortholog_uniprot == row['#ID(s) interactor A'].split(':')[1]),
                            (OrthologEcoli.strain_interactor == strains[num])).one().interactor_id
                        interactors[num].append(session.query(Interactor).filter(Interactor.id == id).one())
                elif (row['#ID(s) interactor A'].split(':')[0] == 'refseq'):
                    if (session.query(OrthologEcoli).filter(
                            (OrthologEcoli.ortholog_refseq == row['#ID(s) interactor A'].split(':')[1]),
                            (OrthologEcoli.strain_interactor == strains[num])).first() != None):
                        id = session.query(OrthologEcoli).filter(
                            (OrthologEcoli.ortholog_refseq == row['#ID(s) interactor A'].split(':')[1]),
                            (OrthologEcoli.strain_interactor == strains[num])).one().interactor_id
                        interactors[num].append(session.query(Interactor).filter(Interactor.id == id).one())
                if (row['ID(s) interactor B'].split(':')[0] == 'uniprotkb'):
                    if (session.query(OrthologEcoli).filter(
                            (OrthologEcoli.ortholog_uniprot == row['ID(s) interactor B'].split(':')[1]),
                            (OrthologEcoli.strain_interactor == strains[num])).first() != None):
                        id = session.query(OrthologEcoli).filter(
                            (OrthologEcoli.ortholog_uniprot == row['ID(s) interactor B'].split(':')[1]),
                            (OrthologEcoli.strain_interactor == strains[num])).one().interactor_id
                        interactors[num].append(session.query(Interactor).filter(Interactor.id == id).one())
                elif (row['ID(s) interactor B'].split(':')[0] == 'refseq'):
                    if (session.query(OrthologEcoli).filter(
                            (OrthologEcoli.ortholog_refseq == row['ID(s) interactor B'].split(':')[1]),
                            (OrthologEcoli.strain_interactor == strains[num])).first() != None):
                        id = session.query(OrthologEcoli).filter(
                            (OrthologEcoli.ortholog_refseq == row['ID(s) interactor B'].split(':')[1]),
                            (OrthologEcoli.strain_interactor == strains[num])).one().interactor_id
                        interactors[num].append(session.query(Interactor).filter(Interactor.id == id).one())

                if (len(interactors[num]) == 2):
                    interaction = session.query(Interaction).filter((Interaction.interactors.contains(interactors[num][0])),
                                                                    (Interaction.interactors.contains(
                                                                        interactors[num][1]))).first()
                    if (interaction != None):
                        interaction.is_Ecoli_ortholog = 2
                    if (interaction == None):
                        interaction = Interaction(strain=strains[num], type=(interactors[num][0].type + '-' + interactors[num][1].type),
                                                  is_experimental = 0, is_Ecoli_ortholog = 1)
                        if (interactors[num][0] == interactors[num][1]):
                            del interactors[num][1]
                        interaction.interactors = interactors[num]
                        session.add(interaction)
                        session.commit()
                    author, date, ref, type, pmid, detection = None, None, None, None, None, None
                    if (row['Publication 1st author(s)'] != '-'):
                        author = row['Publication 1st author(s)'].split('-')[0]
                        date = row['Publication 1st author(s)'].split('-')[1]
                        ref = row['Publication 1st author(s)']
                    if (row['Interaction type(s)'] != '-'):
                        type = row['Interaction type(s)'].split('(')[1][:-1]
                    if (row['Publication Identifier(s)'] != '-'):
                        pmid = row['Publication Identifier(s)'].split('med:')[1][:8]
                    if (row['Interaction detection method(s)'] != '-'):
                        detection = row['Interaction detection method(s)'].split('(')[1][:-1]
                    # there are more than one pmid sometimes
                    reference = InteractionReference(interaction_type=type,
                                                     interaction_id=interaction.id,
                                                     confidence_score=row['Confidence value(s)'],
                                                     pmid=pmid,
                                                     detection_method=detection,
                                                     author_last_name=author,
                                                     publication_date=date,
                                                     publication_ref=ref,
                                                     source_db=row['Source database(s)'].split('(')[1][:-1],
                                                     experimental_role_a=
                                                     row['Experimental role(s) interactor A'].split('(')[1][:-1],
                                                     experimental_role_b=
                                                     row['Experimental role(s) interactor B'].split('(')[1][:-1],
                                                     interactor_a_id=row['#ID(s) interactor A'].split(':')[1],
                                                     interactor_b_id=row['ID(s) interactor B'].split(':')[1])
                    session.add(reference)
        session.commit()
        print(session.query(Interaction).count())

    with open('PSICQUIC_Ecoli_mentha.txt') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
            interactorsPAO1, interactorsPA14 = [], []
            interactors = [interactorsPAO1, interactorsPA14]
            strains = ['PAO1', 'PA14']

            for num in range(0, 2):
                if (session.query(OrthologEcoli).filter(
                        (OrthologEcoli.ortholog_uniprot == row['#ID(s) interactor A'].split(':')[1]),
                        (OrthologEcoli.strain_interactor == strains[num])).first() != None):
                    id = session.query(OrthologEcoli).filter(
                        (OrthologEcoli.ortholog_uniprot == row['#ID(s) interactor A'].split(':')[1]),
                        (OrthologEcoli.strain_interactor == strains[num])).one().interactor_id
                    interactors[num].append(session.query(Interactor).filter(Interactor.id == id).one())
                if (session.query(OrthologEcoli).filter(
                        (OrthologEcoli.ortholog_uniprot == row['ID(s) interactor B'].split(':')[1]),
                        (OrthologEcoli.strain_interactor == strains[num])).first() != None):
                    id = session.query(OrthologEcoli).filter(
                        (OrthologEcoli.ortholog_uniprot == row['ID(s) interactor B'].split(':')[1]),
                        (OrthologEcoli.strain_interactor == strains[num])).one().interactor_id
                    interactors[num].append(session.query(Interactor).filter(Interactor.id == id).one())

                if (len(interactors[num]) == 2):
                    interaction = session.query(Interaction).filter((Interaction.interactors.contains(interactors[num][0])),
                                                                    (Interaction.interactors.contains(
                                                                        interactors[num][1]))).first()
                    if (interaction != None):
                        interaction.is_Ecoli_ortholog = 2
                    if (interaction == None):
                        interaction = Interaction(strain=strains[num], type=(interactors[num][0].type + '-' + interactors[num][1].type),
                                                  is_experimental = 0, is_Ecoli_ortholog = 1)
                        if (interactors[num][0] == interactors[num][1]):
                            del interactors[num][1]
                        interaction.interactors = interactors[num]
                        session.add(interaction)
                        session.commit()
                    reference = InteractionReference(interaction_type=row['Interaction type(s)'].split('(')[1][:-1],
                                                     interaction_id=interaction.id,
                                                     confidence_score=row['Confidence value(s)'],
                                                     pmid=row['Publication Identifier(s)'].split('med:')[1][:8],
                                                     detection_method=
                                                     row['Interaction detection method(s)'].split('(')[1][:-1],
                                                     source_db=row['Source database(s)'].split('(')[1][:-1])
                    session.add(reference)
        session.commit()
        print(session.query(Interaction).count())

    with open('PSICQUIC_Ecoli_MINT.txt') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
            interactorsPAO1, interactorsPA14 = [], []
            interactors = [interactorsPAO1, interactorsPA14]
            strains = ['PAO1', 'PA14']

            for num in range(0, 2):
                if (session.query(OrthologEcoli).filter(
                        (OrthologEcoli.ortholog_uniprot == row['#ID(s) interactor A'].split(':')[1]),
                        (OrthologEcoli.strain_interactor == strains[num])).first() != None):
                    id = session.query(OrthologEcoli).filter(
                        (OrthologEcoli.ortholog_uniprot == row['#ID(s) interactor A'].split(':')[1]),
                        (OrthologEcoli.strain_interactor == strains[num])).one().interactor_id
                    interactors[num].append(session.query(Interactor).filter(Interactor.id == id).one())
                if (session.query(OrthologEcoli).filter(
                        (OrthologEcoli.ortholog_uniprot == row['ID(s) interactor B'].split(':')[1]),
                        (OrthologEcoli.strain_interactor == strains[num])).first() != None):
                    id = session.query(OrthologEcoli).filter(
                        (OrthologEcoli.ortholog_uniprot == row['ID(s) interactor B'].split(':')[1]),
                        (OrthologEcoli.strain_interactor == strains[num])).one().interactor_id
                    interactors[num].append(session.query(Interactor).filter(Interactor.id == id).one())

                if (len(interactors[num]) == 2):
                    interaction = session.query(Interaction).filter(
                        (Interaction.interactors.contains(interactors[num][0])),
                        (Interaction.interactors.contains(
                            interactors[num][1]))).first()
                    if (interaction != None):
                        interaction.is_Ecoli_ortholog = 2
                    if (interaction == None):
                        interaction = Interaction(strain=strains[num],
                                                  type=(interactors[num][0].type + '-' + interactors[num][1].type),
                                                  is_experimental=0, is_Ecoli_ortholog=1)
                        if (interactors[num][0] == interactors[num][1]):
                            del interactors[num][1]
                        interaction.interactors = interactors[num]
                        session.add(interaction)
                        session.commit()
                    reference = InteractionReference(interaction_type=row['Interaction type(s)'].split('(')[1][:-1],
                                                     interaction_id=interaction.id,
                                                     confidence_score=row['Confidence value(s)'],
                                                     pmid=row['Publication Identifier(s)'].split('med:')[1][:8],
                                                     detection_method=row['Interaction detection method(s)'].split('(')[
                                                                          1][
                                                                      :-1],
                                                     author_last_name=row['Publication 1st author(s)'].split(',')[0],
                                                     publication_date=row['Publication 1st author(s)'].split('(')[1][
                                                                      :-1],
                                                     publication_ref=row['Publication 1st author(s)'],
                                                     source_db=row['Source database(s)'].split('(')[1][:-1],
                                                     experimental_role_a=
                                                     row['Experimental role(s) interactor A'].split('(')[1][:-1],
                                                     experimental_role_b=
                                                     row['Experimental role(s) interactor B'].split('(')[1][:-1],
                                                     interactor_a_id=row['#ID(s) interactor A'].split(':')[1],
                                                     interactor_b_id=row['ID(s) interactor B'].split(':')[1])
                    session.add(reference)
        session.commit()
        print(session.query(Interaction).count())

    with open('PSICQUIC_Ecoli_MPIDB.txt') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
            interactorsPAO1, interactorsPA14 = [], []
            interactors = [interactorsPAO1, interactorsPA14]
            strains = ['PAO1', 'PA14']

            for num in range(0, 2):
                if (session.query(OrthologEcoli).filter(
                        (OrthologEcoli.ortholog_uniprot == row['#ID(s) interactor A'].split(':')[1]),
                        (OrthologEcoli.strain_interactor == strains[num])).first() != None):
                    id = session.query(OrthologEcoli).filter(
                        (OrthologEcoli.ortholog_uniprot == row['#ID(s) interactor A'].split(':')[1]),
                        (OrthologEcoli.strain_interactor == strains[num])).one().interactor_id
                    interactors[num].append(session.query(Interactor).filter(Interactor.id == id).one())
                if (session.query(OrthologEcoli).filter(
                        (OrthologEcoli.ortholog_uniprot == row['ID(s) interactor B'].split(':')[1]),
                        (OrthologEcoli.strain_interactor == strains[num])).first() != None):
                    id = session.query(OrthologEcoli).filter(
                        (OrthologEcoli.ortholog_uniprot == row['ID(s) interactor B'].split(':')[1]),
                        (OrthologEcoli.strain_interactor == strains[num])).one().interactor_id
                    interactors[num].append(session.query(Interactor).filter(Interactor.id == id).one())

                if (len(interactors[num]) == 2):
                    interaction = session.query(Interaction).filter(
                        (Interaction.interactors.contains(interactors[num][0])),
                        (Interaction.interactors.contains(
                            interactors[num][1]))).first()
                    if (interaction != None):
                        interaction.is_Ecoli_ortholog = 2
                    if (interaction == None):
                        interaction = Interaction(strain=strains[num],
                                                  type=(interactors[num][0].type + '-' + interactors[num][1].type),
                                                  is_experimental=0, is_Ecoli_ortholog=1)
                        if (interactors[num][0] == interactors[num][1]):
                            del interactors[num][1]
                        interaction.interactors = interactors[num]
                        session.add(interaction)
                        session.commit()
                    reference = InteractionReference(interaction_type=row['Interaction type(s)'].split('(')[1][:-1],
                                                     interaction_id=interaction.id,
                                                     confidence_score=row['Confidence value(s)'],
                                                     pmid=row['Publication Identifier(s)'].split('med:')[1][:8],
                                                     detection_method=row['Interaction detection method(s)'].split('(')[
                                                                          1][
                                                                      :-1],
                                                     author_last_name=row['Publication 1st author(s)'].split(' ')[0],
                                                     publication_date=row['Publication 1st author(s)'].split('(')[1][
                                                                      :-1],
                                                     publication_ref=row['Publication 1st author(s)'],
                                                     source_db=row['Source database(s)'].split('(')[1][:-1],
                                                     experimental_role_a=
                                                     row['Experimental role(s) interactor A'].split('(')[1][:-1],
                                                     experimental_role_b=
                                                     row['Experimental role(s) interactor B'].split('(')[1][:-1],
                                                     interactor_a_id=row['#ID(s) interactor A'].split(':')[1],
                                                     interactor_b_id=row['ID(s) interactor B'].split(':')[1])
                    session.add(reference)
        session.commit()
        print(session.query(Interaction).count())

    with open('DIP_Ecoli.txt') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
            interactorsPAO1, interactorsPA14 = [], []
            interactors = [interactorsPAO1, interactorsPA14]
            strains = ['PAO1', 'PA14']

            for num in range(0, 2):
                ids_A = row['ID Interactor A'].split('|')
                ids_B = row['ID Interactor B'].split('|')
                refseq_A, uniprotkb_A, refseq_B, uniprotkb_B = '', '', '', ''
                for id in ids_A:
                    fields = id.split(':')
                    if (fields[0] == 'refseq'):
                        refseq_A = fields[1]
                    elif (fields[0] == 'uniprotkb'):
                        uniprotkb_A = fields[1]
                for id in ids_B:
                    fields = id.split(':')
                    if (fields[0] == 'refseq'):
                        refseq_B = fields[1]
                    elif (fields[0] == 'uniprotkb'):
                        uniprotkb_B = fields[1]

                has_interactor_A, has_interactor_B = False, False
                if (uniprotkb_A != ''):
                    if (session.query(OrthologEcoli).filter(
                            (OrthologEcoli.ortholog_uniprot == uniprotkb_A),
                            (OrthologEcoli.strain_interactor == strains[num])).first() != None):
                        id = session.query(OrthologEcoli).filter(
                            (OrthologEcoli.ortholog_uniprot == uniprotkb_A),
                            (OrthologEcoli.strain_interactor == strains[num])).one().interactor_id
                        interactors[num].append(session.query(Interactor).filter(Interactor.id == id).one())
                        has_interactor_A = True
                if ((has_interactor_A == False) & (refseq_A != '')):
                    if (session.query(OrthologEcoli).filter(
                            (OrthologEcoli.ortholog_refseq == row['#ID(s) interactor A'].split(':')[1]),
                            (OrthologEcoli.strain_interactor == strains[num])).first() != None):
                        id = session.query(OrthologEcoli).filter(
                            (OrthologEcoli.ortholog_refseq == row['#ID(s) interactor A'].split(':')[1]),
                            (OrthologEcoli.strain_interactor == strains[num])).one().interactor_id
                        interactors[num].append(session.query(Interactor).filter(Interactor.id == id).one())
                if (uniprotkb_B != ''):
                    if (session.query(OrthologEcoli).filter(
                            (OrthologEcoli.ortholog_uniprot == uniprotkb_B),
                            (OrthologEcoli.strain_interactor == strains[num])).first() != None):
                        id = session.query(OrthologEcoli).filter(
                            (OrthologEcoli.ortholog_uniprot == uniprotkb_A),
                            (OrthologEcoli.strain_interactor == strains[num])).one().interactor_id
                        interactors[num].append(session.query(Interactor).filter(Interactor.id == id).one())
                        has_interactor_B = True
                if ((has_interactor_B == False) & (refseq_A != '')):
                    if (session.query(OrthologEcoli).filter(
                            (OrthologEcoli.ortholog_refseq == row['#ID(s) interactor A'].split(':')[1]),
                            (OrthologEcoli.strain_interactor == strains[num])).first() != None):
                        id = session.query(OrthologEcoli).filter(
                            (OrthologEcoli.ortholog_refseq == row['#ID(s) interactor A'].split(':')[1]),
                            (OrthologEcoli.strain_interactor == strains[num])).one().interactor_id
                        interactors[num].append(session.query(Interactor).filter(Interactor.id == id).one())

                if (len(interactors[num]) == 2):
                    interaction = session.query(Interaction).filter((Interaction.interactors.contains(interactors[num][0])),
                                                                    (Interaction.interactors.contains(
                                                                        interactors[num][1]))).first()
                    if (interaction != None):
                        interaction.is_Ecoli_ortholog = 2
                    if (interaction == None):
                        interaction = Interaction(strain=strains[num], type=(interactors[num][0].type + '-' + interactors[num][1].type),
                                                  is_experimental = 0, is_Ecoli_ortholog = 1)
                        if (interactors[num][0] == interactors[num][1]):
                            del interactors[num][1]
                        interaction.interactors = interactors[num]
                        session.add(interaction)
                        session.commit()

                    detections, pmids, types, list = [], [], [], []
                    if (row['Interaction detection method(s)'] != '-'):
                        detections = row['Interaction detection method(s)'].split('|')
                        list.append(detections)
                    if (row['Publication Identifiers'] != '-'):
                        pmids = row['Publication Identifiers'].split('|')
                        list.append(pmids)
                    if (row['Interaction type(s)']):
                        types = row['Interaction type(s)'].split('|')
                        list.append(types)

                    if (len(list) != 0) & (all((len(item) == len(list[0])) for item in list)):
                        for num in range(0, len(list[0])):
                            type, pmid, detection = None, None, None
                            for item in list:
                                if (item == types):
                                    type = types[num].split('(')[1][:-1]
                                if (item == pmids):
                                    pmid = pmids[num*2].split('med:')[1][:8]
                                if (item == detections):
                                    detection = detections[num].split('(')[1][:-1]
                            # there are more than one pmid sometimes
                            reference = InteractionReference(interaction_type=type, interaction_id=interaction.id,
                                                             pmid=pmid, detection_method=detection,
                                                             source_db=row['Source database(s)'].split('(')[1][:-1])
                            session.add(reference)
        session.commit()
        print(session.query(Interaction).count())

    with open('Ecoli_network_tf_gene.txt') as csvfile:
        reader = csv.DictReader(csvfile, delimiter = '\t')
        for row in reader:
            interactorsPAO1, interactorsPA14 = [], []
            interactors = [interactorsPAO1, interactorsPA14]
            strains = ['PAO1', 'PA14']

            for num in range(0, 2):
                if (session.query(OrthologEcoli).filter((OrthologEcoli.ortholog_name == row['TF name'].lower()),
                                                        (OrthologEcoli.strain_interactor == strains[num])).first() != None):
                    id = session.query(OrthologEcoli).filter(
                        (OrthologEcoli.name == row['TF name'].lower()),
                        (OrthologEcoli.strain_interactor == strains[num])).one().interactor_id
                    interactors[num].append(session.query(Interactor).filter(Interactor.id == id).one())
                if (session.query(OrthologEcoli).filter((OrthologEcoli.ortholog_name == row['Regulated gene']),
                                                        (OrthologEcoli.strain_interactor == strains[num])).first() != None):
                    id = session.query(OrthologEcoli).filter(
                        (OrthologEcoli.name == row['Regulated gene'].lower()),
                        (OrthologEcoli.strain_interactor == strains[num])).one().interactor_id
                    interactors[num].append(session.query(Interactor).filter(Interactor.id == id).one())


                if (len(interactors[num]) == 2):
                    # string.append(row['Regulator (TF or sigma)'] + ':' + row['Target'] + '(' + strain + ')')
                    interaction = session.query(Interaction).filter(
                        (Interaction.interactors.contains(interactors[num][0])),
                        (Interaction.interactors.contains(interactors[num][1]))).first()
                    if (interactors[num][0].type != 'CDS(TF/sigma)'):
                        interactors[num][0].type = 'CDS(TF/sigma)'

                    if (interaction != None):
                        if (interaction.type != 'TF/sigma-binding site'):
                            interaction.type = 'TF/sigma-binding site'
                            session.commit()
                        interaction.is_Ecoli_ortholog = 1

                    second_interactor = interactors[num][1]
                    if (interaction == None):
                        interaction = Interaction(strain=strains[num], type='TF/sigma-binding site', is_experimental = 0,
                                                  is_Ecoli_ortholog = 2)
                        if (interactors[0] == interactors[1]):
                            del interactors[1]
                        interaction.interactors = interactors
                        session.add(interaction)
                        session.commit()

                    reference = InteractionReference(interaction_id=interaction.id, pmid=row['pmid'],
                                                     interaction_type='TF/sigma-binding site (' + row['Regulatory effect'] +
                                                                      'regulation)',
                                                     detection_method=row['Evidence'],
                                                     interaction_full_name=interactors[0].id + ' regulates(' +
                                                                           row['Regulatory effect'] + ') ' +
                                                                           second_interactor.id,
                                                     source_db = 'RegulonDB')
                    session.add(reference)
        session.commit()
        print(session.query(InteractionReference).count())