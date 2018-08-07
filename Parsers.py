if __name__ == '__main__':

    import csv
    from Bio.KEGG.REST import kegg_list, kegg_get, kegg_conv
    from Bio.KEGG.KGML.KGML_parser import read
    from Bio.KEGG.KGML.KGML_pathway import Relation
    from Schema1 import Base, Interactor, Metabolite, Protein, ProteinXref, Reference, Localization, GeneOntology, \
        OrthologEcoli, OrthologPseudomonas, Interaction, InteractionReference, InteractionXref
    from sqlalchemy import create_engine, or_
    from sqlalchemy.orm import sessionmaker
    from os.path import exists

    engine = create_engine('sqlite:///:memory:', echo=True)
    Base.metadata.create_all(engine)
    Session = sessionmaker(bind=engine)
    session = Session()
    # dictionary of all compounds in KEGG (kegg_id:[name, pubchem])
    kegg_compounds = {}
    # dictionary of all compounds in ECOCYC (name: [EcoCyc, PubChem, CAS, ChEBI, KEGG])
    ecocyc_compounds = {}

    def find_type(psi_code):
        return {
            '1110': 'predicted interaction',
            '2232': 'molecular association',
            '0914': 'association',
            '0915': 'physical association',
            '0407': 'direct interaction',
        }[psi_code]


    def find_type_KEGG(attrib):
        return {
            'ECrel': 'enzyme-enzyme relation',
            'PPrel': 'protein-protein interaction',
            'GErel': 'gene expression interaction',
            'PCrel': 'protein-compound interaction',
        }[attrib]


    def getGeneInfo(file, strain):
        with open(file) as csvfile:
            fieldnames = ["Sequence", "Locus Tag", "Feature Type", "Start", "End", "Strand", "Name", "Product Name",
                          "Synonyms", "Accession"]
            reader = csv.DictReader(csvfile, fieldnames=fieldnames)
            row_num = 0
            for row in reader:
                if (row['Feature Type'] == 'CDS') & (row_num >= 3):
                    protein = Protein(id=row['Locus Tag'][:-1], name=row['Name'], type='protein', strain=strain,
                                      description=row['Product Name'], accession=row['Accession'])
                    session.add(protein)
                row_num += 1
            session.commit()
        print(session.query(Interactor).count())


    def getXrefs(file, strain):
        with open(file) as csvfile:
            # types = ['RefSeq Accession', 'UniProtKB Accession', 'UniProtKB ID', 'GI Number', 'Uniparc', 'UniRef100 ID',
            #          'UniRef90 ID', 'UniRef50 ID']
            types = ['UniProtKB Accession', 'RefSeq Accession']
            reader = csv.DictReader(csvfile)
            for row in reader:
                if (row['Feature Type'] == 'CDS'):
                    for type in types:
                        if (row[type] == ''): continue
                        if (session.query(ProteinXref).filter((ProteinXref.accession == row[type]),
                                                              (ProteinXref.protein_id == row[
                                                                  'Locus Tag'])).count() == 0):
                            xref = ProteinXref(accession=row[type], protein_id=row['Locus Tag'], data_source=type)
                            session.add(xref)

                        if (type == 'UniProtKB Accession'):
                            if (session.query(Interactor).filter(Interactor.id == row['Locus Tag']).first() != None):
                                interactor = session.query(Interactor).filter(Interactor.id == row['Locus Tag']).one()
                                interactor.uniprotkb = row[type]
                                session.commit()
                            complex_interactors = \
                                session.query(Protein).filter(Protein.uniprotkb == row[type]).count()
                            if (complex_interactors == 2):
                                interactor = Protein(id=row[type], strain=strain, type='protein',
                                                     description='protein complex')
                                session.add(interactor)
                    session.commit()


    def getLocalizations(file):
        with open(file) as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                localization = Localization(protein_id=row['Locus Tag'], localization=row['Subcellular Localization'],
                                            confidence=row['Confidence'])
                session.add(localization)
            session.commit()


    def getOntology(file):
        with open(file) as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                ontology = GeneOntology(protein_id=row['Locus Tag'], accession=row['Accession'], go_term=row['GO Term'],
                                        evidence_code=row['GO Evidence Code'], pmid=row['PMID'],
                                        eco_code=row['Evidence Ontology ECO Code'],
                                        eco_term=row['Evidence Ontology Term'])
                session.add(ontology)
            session.commit()


    def orthologs_ecoli(file, strain):
        # inparalog dictionary, E.coli: P.aeruginosa
        inparalogs = {}

        # get all inparalogs and store in inparalogs dict
        with open(file) as csvfile:
            reader = csv.DictReader(csvfile)

            for row in reader:
                # ignore Non/Borderline SSD
                if (row['Ortholuge Class'] != '') & (row['Ortholuge Class'] != 'SSD'): continue

                if (row['Strain 1 Inparalogs (Locus Tag/Name)'] != ''):
                    strain1_inparalogs = row['Strain 1 Inparalogs (Locus Tag/Name)'].split(']')

                    for inparalog in strain1_inparalogs:
                        if (inparalog == ''): continue
                        trimmed_inparalog = inparalog[:5]
                        if (inparalog[0] == ';'):
                            trimmed_inparalog = inparalog[1:6]
                        if (trimmed_inparalog in inparalogs.keys()):
                            inparalogs[trimmed_inparalog].append(row['Locus Tag (Strain 2)'])
                        else:
                            inparalogs[trimmed_inparalog] = [row['Locus Tag (Strain 2)']]

                if (row['Strain 2 Inparalogs (Locus Tag/Name)'] != ''):
                    strain2_inparalogs = row['Strain 2 Inparalogs (Locus Tag/Name)'].split(']')

                    for inparalog in strain2_inparalogs:
                        if (inparalog == ''): continue
                        trimmed_inparalog = inparalog[:6]
                        if (inparalog[0] == ';'):
                            trimmed_inparalog = inparalog[1:7]
                        if (row['Locus Tag (Strain 1)'] in inparalogs.keys()):
                            inparalogs[row['Locus Tag (Strain 1)']].append(trimmed_inparalog)
                        else:
                            inparalogs[row['Locus Tag (Strain 1)']] = [trimmed_inparalog]

        # parse all orthologs (don't include inparalogs)
        with open(file) as csvfile:
            reader2 = csv.DictReader(csvfile)
            for row in reader2:
                if (row['Ortholuge Class'] != '') & (row['Ortholuge Class'] != 'SSD'): continue

                if (row['Locus Tag (Strain 1)'] in inparalogs.keys()):
                    if (row['Locus Tag (Strain 2)'] in inparalogs[row['Locus Tag (Strain 1)']]):
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

        # get all uniprot ids and gene names for E. coli orthologs
        with open('Ecoli_IDs.csv') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                if (row['Gene names'] != ''):
                    locus = row['Gene names'].split(' ')[0]
                    for ortholog in session.query(OrthologEcoli).filter(OrthologEcoli.ortholog_id == locus).all():
                        ortholog.ortholog_uniprot = row['Entry']
                        ortholog.ortholog_name = row['Gene name']
                    session.commit()

    getGeneInfo('pseudomonas_db_info_PAO1.csv', 'PAO1')
    getGeneInfo('pseudomonas_db_info_PA14.csv', 'PA14')

    getXrefs('pseudomonas_db_xrefs_PAO1.csv', 'PAO1')
    getXrefs('pseudomonas_db_xrefs_PA14.csv', 'PA14')

    #getLocalizations('pseudomonas_db_localizations_PAO1.csv')
    #getLocalizations('pseudomonas_db_localizations_PA14.csv')

    #getOntology('pseudomonas_db_GO_PAO1.csv')
    #getOntology('pseudomonas_db_GO_PA14.csv')

    #orthologs_ecoli('orthologs_PAO1-Ecoli.csv', 'PAO1')
    #orthologs_ecoli('orthologs_PA14-Ecoli.csv', 'PA14')

    def parse_IMEx(file, strain, taxid):
        all_ref =[]
        with open(file) as csvfile:
            reader = csv.DictReader(csvfile, delimiter='\t')
            for row in reader:
                all_ref.append(row['Publication Identifier(s)'])
        #         interactors = []
        #
        #         if ((row['Taxid interactor A'].split('|')[0] != taxid) |
        #                 (row['Taxid interactor B'].split('|')[0] != taxid)): continue
        #
        #         if (session.query(Interactor).filter(
        #                 Interactor.id == row['#ID(s) interactor A'].split(':')[1]).first() != None):
        #             interactors.append(session.query(Interactor).filter(
        #                 Interactor.id == row['#ID(s) interactor A'].split(':')[1]).one())
        #         elif row['#ID(s) interactor A'].split(':')[0] == 'uniprotkb':
        #             if (session.query(Protein).filter(
        #                 Protein.uniprotkb == row['#ID(s) interactor A'].split(':')[1]).first() != None):
        #                 interactors.append(session.query(Protein).filter(
        #                     Protein.uniprotkb == row['#ID(s) interactor A'].split(':')[1]).one())
        #         if (session.query(Interactor).filter(
        #                 Interactor.id == row['ID(s) interactor B'].split(':')[1]).first() != None):
        #             interactors.append(session.query(Interactor).filter(
        #                 Interactor.id == row['ID(s) interactor B'].split(':')[1]).one())
        #         elif row['ID(s) interactor B'].split(':')[0] == 'uniprotkb':
        #             if (session.query(Protein).filter(
        #                 Protein.uniprotkb == row['ID(s) interactor B'].split(':')[1]).first() != None):
        #                 interactors.append(session.query(Protein).filter(
        #                     Protein.uniprotkb == row['ID(s) interactor B'].split(':')[1]).one())
        #
        #         if (len(interactors) != 2): continue
        #
        #         homogenous = (interactors[0] == interactors[1])
        #         interaction = session.query(Interaction).filter((Interaction.interactors.contains(interactors[0])),
        #                                                         (Interaction.interactors.contains(interactors[1])),
        #                                                         (Interaction.homogenous == homogenous)).first()
        #         if (interaction == None):
        #             interaction = Interaction(strain=strain, type=(interactors[0].type + '-' + interactors[1].type),
        #                                       homogenous = homogenous, interactors = interactors)
        #             session.add(interaction), session.commit()
        #
        #         reference = InteractionReference(interaction_type=row['Interaction type(s)'].split('(')[1][:-1],
        #                                          interaction_id=interaction.id,
        #                                          confidence_score=row['Confidence value(s)'],
        #                                          pmid=row['Publication Identifier(s)'].split('med:')[1][:8],
        #                                          detection_method=
        #                                          row['Interaction detection method(s)'].split('(')[1][:-1],
        #                                          author_last_name=row['Publication 1st author(s)'].split(' ')[0],
        #                                          publication_date=row['Publication 1st author(s)'].split('(')[1][:-1],
        #                                          publication_ref=row['Publication 1st author(s)'],
        #                                           # MINT, mpidb, IntAct, DIP
        #                                          source_db=row['Source database(s)'].split('(')[1][:-1],
        #                                          experimental_role_a=
        #                                          row['Experimental role(s) interactor A'].split('(')[1][:-1],
        #                                          experimental_role_b=
        #                                          row['Experimental role(s) interactor B'].split('(')[1][:-1],
        #                                          interactor_a_id=row['#ID(s) interactor A'].split(':')[1],
        #                                          interactor_b_id=row['ID(s) interactor B'].split(':')[1])
        #         interaction.is_experimental = 1
        #         session.add(reference), session.commit()
        # print(session.query(Interaction).count())
        print(all_ref)


    def parse_IntAct(file, strain, taxid):
        with open(file) as csvfile:
            reader = csv.DictReader(csvfile, delimiter='\t')
            for row in reader:
                interactors = []

                if ((row['Taxid interactor A'].split('|')[0] != taxid) |
                        (row['Taxid interactor B'].split('|')[0] != taxid)): continue

                if (session.query(Interactor).filter(
                        Interactor.id == row['#ID(s) interactor A'].split(':')[1]).first() != None):
                    interactors.append(session.query(Interactor).filter(
                        Interactor.id == row['#ID(s) interactor A'].split(':')[1]).one())
                elif row['#ID(s) interactor A'].split(':')[0] == 'uniprotkb':
                    if (session.query(Protein).filter(
                        Protein.uniprotkb == row['#ID(s) interactor A'].split(':')[1]).first() != None):
                        interactors.append(session.query(Protein).filter(
                            Protein.uniprotkb == row['#ID(s) interactor A'].split(':')[1]).one())
                if (session.query(Interactor).filter(
                        Interactor.id == row['ID(s) interactor B'].split(':')[1]).first() != None):
                    interactors.append(session.query(Interactor).filter(
                        Interactor.id == row['ID(s) interactor B'].split(':')[1]).one())
                elif row['ID(s) interactor B'].split(':')[0] == 'uniprotkb':
                    if (session.query(Protein).filter(
                        Protein.uniprotkb == row['ID(s) interactor B'].split(':')[1]).first() != None):
                        interactors.append(session.query(Protein).filter(
                            Protein.uniprotkb == row['ID(s) interactor B'].split(':')[1]).one())

                if (len(interactors) != 2): continue

                homogenous = (interactors[0] == interactors[1])
                interaction = session.query(Interaction).filter((Interaction.interactors.contains(interactors[0])),
                                                                (Interaction.interactors.contains(interactors[1])),
                                                                (Interaction.homogenous == homogenous)).first()
                if (interaction == None):
                    interaction = Interaction(strain=strain, type=(interactors[0].type + '-' + interactors[1].type),
                                              homogenous=homogenous, interactors = interactors)
                    session.add(interaction), session.commit()
                    # how to avoid duplicate references
                reference = InteractionReference(interaction_type=row['Interaction type(s)'].split('(')[1][:-1],
                                                 interaction_id=interaction.id,
                                                 confidence_score=row['Confidence value(s)'],
                                                 pmid=row['Publication Identifier(s)'].split('med:')[1][:8],
                                                 detection_method=
                                                 row['Interaction detection method(s)'].split('(')[1][:-1],
                                                 author_last_name=row['Publication 1st author(s)'].split(',')[0],
                                                 publication_date=row['Publication 1st author(s)'].split('(')[1][:-1],
                                                 publication_ref=row['Publication 1st author(s)'],
                                                 source_db=row['Source database(s)'].split('(')[1][:-1],
                                                 experimental_role_a=
                                                 row['Experimental role(s) interactor A'].split('(')[1][:-1],
                                                 experimental_role_b=
                                                 row['Experimental role(s) interactor B'].split('(')[1][:-1],
                                                 interactor_a_id=row['#ID(s) interactor A'].split(':')[1],
                                                 interactor_b_id=row['ID(s) interactor B'].split(':')[1])
                interaction.is_experimental = 1
                session.add(reference), session.commit()
            print(session.query(Interaction).count())


    def parse_iRefIndex(file, strain, taxid):
        with open(file) as csvfile:
            reader = csv.DictReader(csvfile, delimiter='\t')
            for row in reader:
                interactors = []

                if ((row['Taxid interactor A'].split('|')[0] != taxid) |
                        (row['Taxid interactor B'].split('|')[0] != taxid)): continue

                if (row['#ID(s) interactor A'].split(':')[0] == 'uniprotkb'):
                    if (session.query(Interactor).filter(
                            Interactor.id == row['#ID(s) interactor A'].split(':')[1]).first() != None):
                        interactors.append(session.query(Interactor).filter(
                            Interactor.id == row['#ID(s) interactor A'].split(':')[1]).one())
                    else:
                        if (session.query(Protein).filter(
                                Protein.uniprotkb == row['#ID(s) interactor A'].split(':')[1]).first() != None):
                            interactors.append(session.query(Protein).filter(
                                Protein.uniprotkb == row['#ID(s) interactor A'].split(':')[1]).one())
                elif (row['#ID(s) interactor A'].split(':')[0] == 'refseq'):
                    if (session.query(ProteinXref).filter(ProteinXref.accession ==
                                                             row['#ID(s) interactor A'].split(':')[1]).first() != None):
                        interactors.append(session.query(ProteinXref).filter(ProteinXref.accession ==
                                                                                row['#ID(s) interactor A'].split(
                                                                                    ':')[1]).one().protein)
                if (row['ID(s) interactor B'].split(':')[0] == 'uniprotkb'):
                    if (session.query(Interactor).filter(
                            Interactor.id == row['ID(s) interactor B'].split(':')[1]).first() != None):
                        interactors.append(session.query(Interactor).filter(
                            Interactor.id == row['ID(s) interactor B'].split(':')[1]).one())
                    elif (session.query(Protein).filter(
                            Protein.uniprotkb == row['ID(s) interactor B'].split(':')[1]).first() != None):
                        interactors.append(session.query(Protein).filter(
                            Protein.uniprotkb == row['ID(s) interactor B'].split(':')[1]).one())
                elif (row['ID(s) interactor B'].split(':')[0] == 'refseq'):
                    if (session.query(ProteinXref).filter(ProteinXref.accession ==
                                                             row['ID(s) interactor B'].split(':')[1]).first() != None):
                        interactors.append(session.query(ProteinXref).filter(ProteinXref.accession ==
                                                                                row['ID(s) interactor B'].split(
                                                                                    ':')[1]).one().protein)

                if (len(interactors) != 2): continue

                homogenous = (interactors[0] == interactors[1])
                interaction = session.query(Interaction).filter((Interaction.interactors.contains(interactors[0])),
                                                                (Interaction.interactors.contains(interactors[1])),
                                                                (Interaction.homogenous == homogenous)).first()
                if (interaction == None):
                    interaction = Interaction(strain=strain, type=(interactors[0].type + '-' + interactors[1].type),
                                              homogenous=homogenous, interactors = interactors)
                    session.add(interaction), session.commit()

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
                reference = InteractionReference(interaction_type=type, interaction_id=interaction.id, pmid=pmid,
                                                 confidence_score=row['Confidence value(s)'], author_last_name=author,
                                                 detection_method=detection, publication_date=date, publication_ref=ref,
                                                 source_db=row['Source database(s)'].split('(')[1][:-1],
                                                 experimental_role_a=
                                                 row['Experimental role(s) interactor A'].split('(')[1][:-1],
                                                 experimental_role_b=
                                                 row['Experimental role(s) interactor B'].split('(')[1][:-1],
                                                 interactor_a_id=row['#ID(s) interactor A'].split(':')[1],
                                                 interactor_b_id=row['ID(s) interactor B'].split(':')[1])
                interaction.is_experimental = 1
                session.add(reference), session.commit()
            print(session.query(Interaction).count())


    def parse_mentha(file, strain, taxid):
        with open(file) as csvfile:
            reader = csv.DictReader(csvfile, delimiter='\t')
            for row in reader:
                interactors = []

                if ((row['Taxid interactor A'].split('|')[0] != taxid) |
                        (row['Taxid interactor B'].split('|')[0] != taxid)): continue

                if (session.query(Interactor).filter(
                        Interactor.id == row['#ID(s) interactor A'].split(':')[1]).first() != None):
                    interactors.append(session.query(Interactor).filter(
                        Interactor.id == row['#ID(s) interactor A'].split(':')[1]).one())
                elif row['#ID(s) interactor A'].split(':')[0] == 'uniprotkb':
                    if (session.query(Protein).filter(
                        Protein.uniprotkb == row['#ID(s) interactor A'].split(':')[1]).first() != None):
                        interactors.append(session.query(Protein).filter(
                            Protein.uniprotkb == row['#ID(s) interactor A'].split(':')[1]).one())
                if (session.query(Interactor).filter(
                        Interactor.id == row['ID(s) interactor B'].split(':')[1]).first() != None):
                    interactors.append(session.query(Interactor).filter(
                        Interactor.id == row['ID(s) interactor B'].split(':')[1]).one())
                elif row['ID(s) interactor B'].split(':')[0] == 'uniprotkb':
                    if (session.query(Protein).filter(
                        Protein.uniprotkb == row['ID(s) interactor B'].split(':')[1]).first() != None):
                        interactors.append(session.query(Protein).filter(
                            Protein.uniprotkb == row['ID(s) interactor B'].split(':')[1]).one())

                if (len(interactors) != 2): continue

                homogenous = (interactors[0] == interactors[1])
                interaction = session.query(Interaction).filter((Interaction.interactors.contains(interactors[0])),
                                                                (Interaction.interactors.contains(interactors[1])),
                                                                (Interaction.homogenous == homogenous)).first()
                if (interaction == None):
                    interaction = Interaction(strain=strain, type=(interactors[0].type + '-' + interactors[1].type),
                                              homogenous=homogenous, interactors = interactors)
                    interaction.interactors = interactors
                    session.add(interaction), session.commit()

                reference = InteractionReference(interaction_type=row['Interaction type(s)'].split('(')[1][:-1],
                                                 interaction_id=interaction.id,
                                                 confidence_score=row['Confidence value(s)'],
                                                 pmid=row['Publication Identifier(s)'].split('med:')[1][:8],
                                                 detection_method=
                                                 row['Interaction detection method(s)'].split('(')[1][:-1],
                                                 source_db=row['Source database(s)'].split('(')[1][:-1])
                interaction.is_experimental = 1
                session.add(reference), session.commit()
            print(session.query(Interaction).count())


    def parse_MINT(file, strain, taxid):
        with open(file) as csvfile:
            reader = csv.DictReader(csvfile, delimiter='\t')
            for row in reader:
                interactors = []

                if ((row['Taxid interactor A'].split('|')[0] != taxid) |
                        (row['Taxid interactor B'].split('|')[0] != taxid)): continue

                if (session.query(Interactor).filter(
                        Interactor.id == row['#ID(s) interactor A'].split(':')[1]).first() != None):
                    interactors.append(session.query(Interactor).filter(
                        Interactor.id == row['#ID(s) interactor A'].split(':')[1]).one())
                elif row['#ID(s) interactor A'].split(':')[0] == 'uniprotkb':
                    if (session.query(Protein).filter(
                            Protein.uniprotkb == row['#ID(s) interactor A'].split(':')[1]).first() != None):
                        interactors.append(session.query(Protein).filter(
                            Protein.uniprotkb == row['#ID(s) interactor A'].split(':')[1]).one())
                if (session.query(Interactor).filter(
                        Interactor.id == row['ID(s) interactor B'].split(':')[1]).first() != None):
                    interactors.append(session.query(Interactor).filter(
                        Interactor.id == row['ID(s) interactor B'].split(':')[1]).one())
                elif row['ID(s) interactor B'].split(':')[0] == 'uniprotkb':
                    if (session.query(Protein).filter(
                            Protein.uniprotkb == row['ID(s) interactor B'].split(':')[1]).first() != None):
                        interactors.append(session.query(Protein).filter(
                            Protein.uniprotkb == row['ID(s) interactor B'].split(':')[1]).one())

                if (len(interactors) != 2): continue

                homogenous = (interactors[0] == interactors[1])
                interaction = session.query(Interaction).filter((Interaction.interactors.contains(interactors[0])),
                                                                (Interaction.interactors.contains(interactors[1])),
                                                                (Interaction.homogenous == homogenous)).first()
                if (interaction == None):
                    interaction = Interaction(strain=strain, type=(interactors[0].type + '-' + interactors[1].type),
                                              homogenous=homogenous, interactors = interactors)
                    session.add(interaction), session.commit()

                reference = InteractionReference(interaction_type=row['Interaction type(s)'].split('(')[1][:-1],
                                                 interaction_id=interaction.id,
                                                 confidence_score=row['Confidence value(s)'],
                                                 pmid=row['Publication Identifier(s)'].split('med:')[1][:8],
                                                 detection_method=
                                                 row['Interaction detection method(s)'].split('(')[1][:-1],
                                                 author_last_name=row['Publication 1st author(s)'].split(',')[0],
                                                 publication_date=row['Publication 1st author(s)'].split('(')[1][:-1],
                                                 publication_ref=row['Publication 1st author(s)'],
                                                 source_db=row['Source database(s)'].split('(')[1][:-1],
                                                 experimental_role_a=
                                                 row['Experimental role(s) interactor A'].split('(')[1][:-1],
                                                 experimental_role_b=
                                                 row['Experimental role(s) interactor B'].split('(')[1][:-1],
                                                 interactor_a_id=row['#ID(s) interactor A'].split(':')[1],
                                                 interactor_b_id=row['ID(s) interactor B'].split(':')[1])
                interaction.is_experimental = 1
                session.add(reference)
            session.commit()
            print(session.query(Interaction).count())


    def get_kegg_compounds():
        # fill in compounds dictionary (kegg_id: [name, pubchem])
        for compound in kegg_list(database='compound'):
            kegg_compounds[compound[4:10]] = {}
            kegg_compounds[compound[4:10]]["name"] = compound.split('\t')[1].split(';')[0]
        for cpd_id in kegg_conv('pubchem', 'compound').read().split('cpd:'):
            if cpd_id != '':
                kegg_compounds[cpd_id[:6]]["PubChem"] = cpd_id.split('pubchem:')[1][:-1]


    def parse_kegg(org_id, strain):
        get_kegg_compounds()
        # get pathways for organism specified by org_id
        pathways = kegg_list(database='pathway', org=org_id).read().split('path:')
        path_names, path_ids = [], []

        for path in pathways:
            if path != '':
                path_names.append(path.split('\t')[1].split(' -')[0])
                path_ids.append(path[:8])

        mas_compound_ids = []
        for path in path_ids:
            # get kgml representation of path
            kgml_path = read(kegg_get(path, option='kgml'))
            path_name = kgml_path._getname()
            # dictionary of compounds in current path (node_id: kegg_id)
            #   compound._getid() returns node id (only relevant in context of current path)
            #   compound._getname() returns kegg id (relevant in overall KEGG DB
            compound_ids = {}
            for compound in kgml_path.compounds:
                compound_ids[compound._getid()] = compound._getname()[-6:]
            mas_compound_ids.append(compound_ids)

            # go through each relation in path
            for relation in kgml_path.relations:
                relation_type = relation.element.attrib['type']

                if (relation_type == 'maplink'): continue
                # relation._getentry1/2() returns  protein id (locus) or compound id (KEGG id)
                entries = [relation._getentry1()._getname(), relation._getentry2()._getname()]
                if (entries[0] == 'undefined') | (entries[1] == 'undefined'): continue
                interactors = [[], []]
                new_metabolites = [[], []]
                # go through each entry in the relation, find/create interactors
                for num in range(0, 2):
                    # each entry may contain >1 id; go through all of them
                    for id in entries[num].split(' '):
                        if (id == ''): continue
                        # if interactor is not protein or compound, continue
                        if (id.split(':')[0] != org_id) & (id.split(':')[1] not in kegg_compounds): continue

                        # check if interactor (protein or metabolite) already exists
                        if (session.query(Interactor).filter(Interactor.id == id.split(':')[1]).first() is not None):
                            interactors[num].append(
                                session.query(Interactor).filter(Interactor.id == id.split(':')[1]).one())
                        # if it doesnt exist, it's not a valid protein, so check if it is a valid compound
                        elif (id.split(':')[1] in kegg_compounds):
                            # if it is a valid compound, create new metabolite
                            new_metabolites[num].append(id.split(':')[1])
                        # if parsing E. coli path, add all orthologs to interactor list
                        elif (org_id == 'eco'):
                            for ortholog in (session.query(OrthologEcoli).filter(
                                    (OrthologEcoli.ortholog_id == id.split(':')[1]),
                                    (OrthologEcoli.strain_protein == strain)).all()):
                                interactors[num].append(ortholog.protein)

                # create list of interactor pairs from two separate lists (interactors[0], interactors[1])
                interactor_pairs = []
                for interactor1 in interactors[0]:
                    for interactor2 in interactors[1]:
                        if (interactor1.type != 'metabolite') | (interactor2.type != 'metabolite'):
                            interactor_pairs.append([interactor1, interactor2])

                for interactor1 in interactors[0]:
                    for interactor2 in new_metabolites[1]:
                        if (interactor1.type != 'metabolite'):
                            new_metabolite = Metabolite(id=interactor2, type='metabolite',
                                                        name=kegg_compounds[interactor2]["name"],
                                                        PubChem=kegg_compounds[interactor2]["PubChem"],
                                                        KEGG=interactor2)
                            session.add(new_metabolite), session.commit()
                            interactor_pairs.append([interactor1, new_metabolite])

                for interactor1 in interactors[1]:
                    for interactor2 in new_metabolites[0]:
                        if (interactor1.type != 'metabolite'):
                            new_metabolite = session.query(Metabolite).filter(Metabolite.id == interactor2).first()
                            if (new_metabolite is None):
                                new_metabolite = Metabolite(id=interactor2, type='metabolite',
                                                        name=kegg_compounds[interactor2]["name"],
                                                        PubChem=kegg_compounds[interactor2]["PubChem"],
                                                        KEGG=interactor2)
                                session.add(new_metabolite), session.commit()
                            interactor_pairs.append([new_metabolite, interactor1])

                # where should this be (include interactions where 1 of primary interactors has no ortholog?)
                if (len(interactor_pairs) == 0): continue

                # get all intermediates in reaction of type compound
                intermeds = []
                for subtype in relation.element.iter(tag='subtype'):
                    if 'compound' in subtype.attrib:
                        compound_node_id = subtype.attrib['compound']
                        if (compound_node_id == None): continue
                        if (int(compound_node_id) not in compound_ids): continue
                        # if compound id is valid, either add existing matching metabolite or create new one and add
                        compound_id = compound_ids[int(compound_node_id)]
                        metabolite = session.query(Metabolite).filter(Metabolite.KEGG == compound_id).first()
                        if metabolite is None:
                            metabolite = Metabolite(KEGG=compound_id, name=kegg_compounds[compound_id]["name"],
                                                    PubChem=kegg_compounds[compound_id]["PubChem"],
                                                    type='metabolite', id=compound_id)
                            session.add(metabolite), session.commit()
                        else:
                            metabolite = session.query(Metabolite).filter(Metabolite.KEGG == compound_id).one()
                        intermeds.append(metabolite)

                # add protein - intermediate interactor pairs
                for interactor_list in interactors:
                    for interactor in interactor_list:
                        if (interactor.type != 'metabolite'):
                            for intermed in intermeds:
                                interactor_pairs.append([interactor, intermed])

                for interactor_pair in interactor_pairs:
                    homogenous = (interactor_pair[0] == interactor_pair[1])
                    interaction = session.query(Interaction).filter(
                        (Interaction.interactors.contains(interactor_pair[0])),
                        (Interaction.interactors.contains(interactor_pair[1])),
                        (Interaction.homogenous == homogenous)).first()

                    if (interaction == None):
                        interaction = Interaction(type=interactor_pair[0].type + '-' + interactor_pair[1].type,
                                                  strain=strain, is_experimental=0, homogenous=homogenous,
                                                  interactors=interactor_pair)
                        if org_id == 'eco':
                            interaction.ortholog_derived = 'from E. coli'
                        session.add(interaction), session.commit()
                    else:
                        interaction = session.query(Interaction).filter(
                            (Interaction.interactors.contains(interactor_pair[0])),
                            (Interaction.interactors.contains(interactor_pair[1])),
                            (Interaction.homogenous == homogenous)).one()
                        if org_id == 'eco':
                            if (interaction.ortholog_derived == None):
                                interaction.ortholog_derived = 'confirmed from E. coli'
                            else:
                                if ('from E.coli' not in interaction.ortholog_derived):
                                    interaction.ortholog_derived += ', confirmed from E.coli'

                    if (session.query(InteractionXref).filter((InteractionXref.interaction_id == interaction.id),
                                                              (InteractionXref.data_source == 'KEGG')).first() == None):
                        xref = InteractionXref(interaction_id=interaction.id, data_source='KEGG')
                        session.add(xref), session.commit()

        print(session.query(Interaction).count())


    def parse_pseudomonas():
        with open('PAO1_GeoffWinsor.csv') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                interactors = []
                if (session.query(Interactor).filter(Interactor.id == row['locus_tag']).first() != None):
                    interactors.append(session.query(Interactor).filter(Interactor.id == row['locus_tag']).one())
                row = next(reader)
                if (session.query(Interactor).filter(Interactor.id == row['locus_tag']).first() != None):
                    interactors.append(session.query(Interactor).filter(Interactor.id == row['locus_tag']).one())

                if (len(interactors) != 2): continue
                homogenous = (interactors[0] == interactors[1])
                interaction = session.query(Interaction).filter((Interaction.interactors.contains(interactors[0])),
                                                                (Interaction.interactors.contains(interactors[1])),
                                                                (Interaction.homogenous == homogenous)).first()
                if (interaction == None):
                    interaction = Interaction(strain='PAO1', type=(interactors[0].type + '-' + interactors[1].type),
                                              homogenous=homogenous, interactors = interactors)
                    session.add(interaction), session.commit()

                reference = InteractionReference(interaction_type=row['type'], interaction_id=interaction.id,
                                                 pmid=row['pmid'], interaction_full_name=row['full_name'],
                                                 detection_method=row['experimental_type'])
                interaction.is_experimental = 1
                session.add(reference)
            session.commit()
            print(session.query(Interaction).count())

        # with open('PAO1_STRING.txt') as csvfile:
        #     fieldnames = ['interactor_A', 'interactor_B', 'altID_A', 'altID_B', 'alias_A', 'alias_B', 'detection',
        #                   'publication', 'publication_ID', 'taxid_A', 'taxid_B', 'type', 'source_db', 'identifier',
        #                   'confidence']
        #     reader = csv.DictReader(csvfile, delimiter='\t', fieldnames=fieldnames)
        #     for row in reader:
        #         interactors = []
        #         locus_tag1 = row['interactor_A'].split('|')[0].split('.')[1]
        #         if (session.query(Interactor).filter(Interactor.id == locus_tag1).first() != None):
        #             interactors.append(session.query(Interactor).filter(Interactor.id == locus_tag1).one())
        #         locus_tag2 = row['interactor_B'].split('|')[0].split('.')[1]
        #         if (session.query(Interactor).filter(Interactor.id == locus_tag2).first() != None):
        #             interactors.append(session.query(Interactor).filter(Interactor.id == locus_tag2).one())
        #
        #         if (len(interactors) != 2): continue
        #         homogenous = (interactors[0] == interactors[1])
        #         interaction = session.query(Interaction).filter((Interaction.interactors.contains(interactors[0])),
        #                                                         (Interaction.interactors.contains(interactors[1])),
        #                                                         (Interaction.homogenous == homogenous)).first()
        #         if (interaction == None):
        #             interaction = Interaction(strain='PAO1', type=(interactors[0].type + '-' + interactors[1].type),
        #                                       homogenous=homogenous, interactors = interactors)
        #             session.add(interaction), session.commit()
        #
        #         detection = row['detection'].split('(')[1][:-1]
        #         confidence = row['confidence'].split(':')[1]
        #         type = find_type(row['type'].split(':')[2][:-1])
        #         source_db, author, date, ref, pmid = None, None, None, None, None
        #
        #         if (row['source_db'] != '-'):
        #             source_db = row['source_db'].split('(')[1][:-1]
        #         if (row['publication'] != '-'):
        #             author = row['publication'].split(' ')[0]
        #             date = row['publication'].split('(')[1][:-1]
        #             ref = row['publication']
        #         if (row['publication_ID'] != '-'):
        #             pmid = row['publication_ID'].split(':')[1]
        #         reference = InteractionReference(interaction_type=type, interaction_id=interaction.id,
        #                                          detection_method=detection, confidence_score=confidence,
        #                                          source_db=source_db, pmid=pmid, author_last_name=author,
        #                                          publication_ref=ref, publication_date=date)
        #         session.add(reference)
        #         if (detection == 'experimental interaction detection'):
        #             interaction.is_experimental = 1
        #     session.commit()
        #     print(session.query(Interaction).count())

        with open('PAO1_Interactome_predicted.csv') as csvfile:
            reader = csv.DictReader(csvfile)

            for row in reader:
                if (float(row['Confidence']) < 0.9): continue
                interactors = []
                if (session.query(Interactor).filter(Interactor.id == row['Protein1']).first() != None):
                    interactors.append(session.query(Interactor).filter(Interactor.id == row['Protein1']).one())

                if (session.query(Interactor).filter(Interactor.id == row['Protein2']).first() != None):
                    interactors.append(session.query(Interactor).filter(Interactor.id == row['Protein2']).one())

                if (len(interactors) != 2): continue
                homogenous = (interactors[0] == interactors[1])
                interaction = session.query(Interaction).filter((Interaction.interactors.contains(interactors[0])),
                                                                (Interaction.interactors.contains(interactors[1])),
                                                                (Interaction.homogenous == homogenous)).first()
                if (interaction == None):
                    interaction = Interaction(strain='PAO1', type=(interactors[0].type + '-' + interactors[1].type),
                                              homogenous=homogenous, is_experimental = 0, interactors = interactors)
                    session.add(interaction), session.commit()

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
                elif (session.query(Protein).filter(Protein.uniprotkb == row['proA']).first() != None):
                    interactors.append(session.query(Protein).filter(Protein.uniprotkb == row['proA']).one())
                if (session.query(Interactor).filter(Interactor.id == row['proB']).first() != None):
                    interactors.append(session.query(Interactor).filter(Interactor.id == row['proB']).one())
                elif (session.query(Protein).filter(Protein.uniprotkb == row['proB']).first() != None):
                    interactors.append(session.query(Protein).filter(Protein.uniprotkb == row['proB']).one())

                if (len(interactors) != 2): continue
                homogenous = (interactors[0] == interactors[1])
                interaction = session.query(Interaction).filter((Interaction.interactors.contains(interactors[0])),
                                                                (Interaction.interactors.contains(interactors[1])),
                                                                (Interaction.homogenous == homogenous)).first()
                if (interaction == None):
                    interaction = Interaction(strain='PAO1', type=(interactors[0].type + '-' + interactors[1].type),
                                              homogenous=homogenous, interactors = interactors)
                    session.add(interaction), session.commit()
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

                if ((row['Taxid interactor A'] != 'taxid:208964(Pseudomonas aeruginosa PAO1)') |
                        (row['Taxid interactor B'] != 'taxid:208964(Pseudomonas aeruginosa PAO1)')): continue

                if (session.query(Interactor).filter(
                        Interactor.id == row['#ID(s) interactor A'].split(':')[1]).first() != None):
                    interactors.append(session.query(Interactor).filter(
                        Interactor.id == row['#ID(s) interactor A'].split(':')[1]).one())
                elif (session.query(Protein).filter(
                        Protein.uniprotkb == row['#ID(s) interactor A'].split(':')[1]).first() != None):
                    interactors.append(session.query(Protein).filter(
                        Protein.uniprotkb == row['#ID(s) interactor A'].split(':')[1]).one())
                if (session.query(Interactor).filter(
                        Interactor.id == row['ID(s) interactor B'].split(':')[1]).first() != None):
                    interactors.append(session.query(Interactor).filter(
                        Interactor.id == row['ID(s) interactor B'].split(':')[1]).one())
                elif (session.query(Protein).filter(
                        Protein.uniprotkb == row['ID(s) interactor B'].split(':')[1]).first() != None):
                    interactors.append(session.query(Protein).filter(
                        Protein.uniprotkb == row['ID(s) interactor B'].split(':')[1]).one())

                if (len(interactors) != 2): continue
                homogenous = (interactors[0] == interactors[1])
                interaction = session.query(Interaction).filter((Interaction.interactors.contains(interactors[0])),
                                                                (Interaction.interactors.contains(interactors[1])),
                                                                (Interaction.homogenous == homogenous)).first()
                if (interaction == None):
                    interaction = Interaction(strain='PAO1', type=(interactors[0].type + '-' + interactors[1].type),
                                              homogenous=homogenous, interactors = interactors)
                    session.add(interaction), session.commit()

                reference = InteractionReference(interaction_type=row['Interaction type(s)'],
                                                 interaction_id=interaction.id,
                                                 confidence_score=row['Confidence value(s)'],
                                                 pmid=row['Publication Identifier(s)'].split('med:')[1][:8],
                                                 detection_method=
                                                 row['Interaction detection method(s)'].split('(')[1][:-1],
                                                 author_last_name=row['Publication 1st author(s)'].split(',')[0],
                                                 publication_date=
                                                 row['Publication 1st author(s)'].split('(')[1][:-1],
                                                 publication_ref=row['Publication 1st author(s)'])
                session.add(reference)
                interaction.is_experimental = 1
                # xref = session.query(InteractionXref).filter(
                #     InteractionXref.accession == row['Interaction identifier(s)'].split(':')[1]).first()
                # if (xref == None):
                #     xref = InteractionXref(accession=row['Interaction identifier(s)'].split(':')[1],
                #                            data_source=row['Source database(s)'].split('(')[1][:-1],
                #                            interaction_id=interaction.id)
                #     session.add(xref)
            session.commit()
            print(session.query(Interaction).count())

        with open('PSICQUIC_Paeruginosa_PAO1_MPIDB.txt') as csvfile:
            reader = csv.DictReader(csvfile, delimiter='\t')
            for row in reader:
                interactors = []

                if ((row['Taxid interactor A'].split('|')[0] != 'taxid:208964(pseae)') |
                        (row['Taxid interactor B'].split('|')[0] != 'taxid:208964(pseae)')): continue

                if (session.query(Interactor).filter(
                        Interactor.id == row['#ID(s) interactor A'].split(':')[1]).first() != None):
                    interactors.append(session.query(Interactor).filter(
                        Interactor.id == row['#ID(s) interactor A'].split(':')[1]).one())
                elif row['#ID(s) interactor A'].split(':')[0] == 'uniprotkb':
                    if (session.query(Protein).filter(
                            Protein.uniprotkb == row['#ID(s) interactor A'].split(':')[1]).first() != None):
                        interactors.append(session.query(Protein).filter(
                            Protein.uniprotkb == row['#ID(s) interactor A'].split(':')[1]).one())
                if (session.query(Interactor).filter(
                        Interactor.id == row['ID(s) interactor B'].split(':')[1]).first() != None):
                    interactors.append(session.query(Interactor).filter(
                        Interactor.id == row['ID(s) interactor B'].split(':')[1]).one())
                elif row['ID(s) interactor B'].split(':')[0] == 'uniprotkb':
                    if (session.query(Protein).filter(
                            Protein.uniprotkb == row['ID(s) interactor B'].split(':')[1]).first() != None):
                        interactors.append(session.query(Protein).filter(
                            Protein.uniprotkb == row['ID(s) interactor B'].split(':')[1]).one())

                if (len(interactors) != 2): continue

                homogenous = (interactors[0] == interactors[1])
                interaction = session.query(Interaction).filter((Interaction.interactors.contains(interactors[0])),
                                                                (Interaction.interactors.contains(interactors[1])),
                                                                (Interaction.homogenous == homogenous)).first()
                if (interaction == None):
                    interaction = Interaction(strain='PAO1', type=(interactors[0].type + '-' + interactors[1].type),
                                              homogenous=homogenous, interactors = interactors)
                    session.add(interaction), session.commit()

                reference = InteractionReference(interaction_type=row['Interaction type(s)'].split('(')[1][:-1],
                                                 interaction_id=interaction.id,
                                                 confidence_score=row['Confidence value(s)'],
                                                 pmid=row['Publication Identifier(s)'].split('med:')[1][:8],
                                                 detection_method=
                                                 row['Interaction detection method(s)'].split('(')[1][:-1],
                                                 author_last_name=row['Publication 1st author(s)'].split(' ')[0],
                                                 publication_date=row['Publication 1st author(s)'].split('(')[1][:-1],
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

        with open('PAO1-PA14_regulatory_network.csv') as csvfile:
            reader = csv.DictReader(csvfile)
            string = []
            for row in reader:
                strains = row['Strain'].split(',')
                for strain in strains:
                    if ((strain != 'PAO1') & (strain != 'PA14')): continue
                    interactors = []
                    if (session.query(Protein).filter((Protein.name == row['Regulator (TF or sigma)']),
                                                         (Protein.strain == strain)).first() is not None):
                        interactors.append(
                            session.query(Protein).filter((Protein.name == row['Regulator (TF or sigma)']),
                                                             (Protein.strain == strain)).one())
                    elif (session.query(Interactor).filter(
                            Interactor.id == row['Regulator (TF or sigma)']).first() is not None):
                            interactors.append(session.query(Interactor).filter(
                                    Interactor.id == row['Regulator (TF or sigma)']).one())
                    if (session.query(Protein).filter((Protein.name == row['Target']),
                                                         (Protein.strain == strain)).first() is not None):
                        interactors.append(session.query(Protein).filter((Protein.name == row['Target']),
                                                                            (Protein.strain == strain)).one())
                    elif (session.query(Interactor).filter(Interactor.id == row['Target']).first() is not None):
                            interactors.append(
                                session.query(Interactor).filter(Interactor.id == row['Target']).one())

                    string.append(len(interactors))
                    if (len(interactors) != 2): continue

                    # string.append(row['Regulator (TF or sigma)'] + ':' + row['Target'] + '(' + strain + ')')
                    homogenous = (interactors[0] == interactors[1])
                    interaction = session.query(Interaction).filter((Interaction.interactors.contains(interactors[0])),
                                                                    (Interaction.interactors.contains(interactors[1])),
                                                                    (Interaction.homogenous == homogenous)).first()
                    interactors[0].is_TF = 1

                    if (interaction == None):
                        interaction = Interaction(strain=strain, type=interactors[0].type + ' - ' + interactors[1].type,
                                                  homogenous = homogenous, interactors = interactors)
                        session.add(interaction), session.commit()

                    reference = InteractionReference(interaction_id=interaction.id, pmid=row['pmid'],
                                                     interaction_type='TF/sigma-binding site (' + row['mode'] +
                                                                      'regulation)',
                                                     detection_method=row['evidence'], source_db=row['source_db'],
                                                     interaction_full_name=interactors[0].id + ' regulates(' +
                                                                           row['mode'] + ') ' + interactors[1].id)
                    session.add(reference)
                    interaction.is_experimental = 1
            session.commit()
            print(session.query(InteractionReference).count())

    parse_IMEx('PSICQUIC_Paeruginosa_PAO1_IMEx.txt', 'PAO1', 'taxid:208964(pseae)')
    parse_IMEx('PSICQUIC_Paeruginosa_PA14_IMEx.txt', 'PA14', 'taxid:208963(pseab)')

    parse_iRefIndex('PSICQUIC_Paeruginosa_PAO1_iRefIndex.txt', 'PAO1', 'taxid:208964(Pseudomonas aeruginosa PAO1)')
    parse_iRefIndex('PSICQUIC_Paeruginosa_PA14_iRefIndex.txt', 'PA14', 'taxid:208963(Pseudomonas aeruginosa UCBPP-PA14)')

    parse_mentha('PSICQUIC_Paeruginosa_PAO1_mentha.txt', 'PAO1', 'taxid:208964(Pseudomonas aeruginosa PAO1)')
    parse_mentha('PSICQUIC_Paeruginosa_PA14_mentha.txt', 'PA14', 'taxid:208963(Pseudomonas aeruginosa UCBPP-PA14)')

    parse_MINT('PSICQUIC_Paeruginosa_PAO1_MINT.txt', 'PAO1', 'taxid:208964(pseae)')
    parse_MINT('PSICQUIC_Paeruginosa_PA14_MINT.txt', 'PA14', 'taxid:208963(pseab)')

    parse_IntAct('PSICQUIC_Paeruginosa_PAO1_IntAct.txt', 'PAO1', 'taxid:208964(pseae)')
    parse_IntAct('PA14_IntAct.txt', 'PA14', 'taxid:208963(pseab)')
    parse_IntAct('PAO1_IntAct.txt', 'PAO1', 'taxid:208964(pseae)')

    parse_kegg('pae', 'PAO1')
    parse_kegg('pau', 'PA14')

    parse_pseudomonas()

    def parse_pseudomonas_orthologs():
        inparalogs = {}
        with open('orthologs_PAO1-PA14.csv') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                if (row['Strain 1 Inparalogs (Locus Tag/Name)'] != ''):
                    strain1_inparalogs = row['Strain 1 Inparalogs (Locus Tag/Name)'].split(']')

                    for inparalog in strain1_inparalogs:
                        if inparalog == '': continue
                        trimmed_inparalog = inparalog[:6]
                        if inparalog[0] == ';':
                            trimmed_inparalog = inparalog[1:7]

                        if trimmed_inparalog in inparalogs.keys():
                            inparalogs[trimmed_inparalog].append(row['Locus Tag (Strain 2)'])
                        else:
                            inparalogs[trimmed_inparalog] = [row['Locus Tag (Strain 2)']]

                if (row['Strain 2 Inparalogs (Locus Tag/Name)'] != ''):
                    strain2_inparalogs = row['Strain 2 Inparalogs (Locus Tag/Name)'].split(']')

                    for inparalog in strain2_inparalogs:
                        if inparalog == '': continue
                        trimmed_inparalog = inparalog[:9]
                        if inparalog[0] == ';':
                            trimmed_inparalog = inparalog[1:10]

                        if row['Locus Tag (Strain 1)'] in inparalogs.keys():
                            inparalogs[row['Locus Tag (Strain 1)']].append(trimmed_inparalog)
                        else:
                            inparalogs[row['Locus Tag (Strain 1)']] = [trimmed_inparalog]

        with open('orthologs_PAO1-PA14.csv') as csvfile:
            reader2 = csv.DictReader(csvfile)
            for row in reader2:
                if (row['Locus Tag (Strain 1)']) in inparalogs.keys():
                    if (row['Locus Tag (Strain 2)']) in inparalogs[row['Locus Tag (Strain 1)']]:
                        continue
                if (session.query(Interactor).filter(Interactor.id == row['Locus Tag (Strain 1)']).first() != None):
                    if (session.query(Interactor).filter(Interactor.id == row['Locus Tag (Strain 2)']).first() != None):
                        ortholog1 = OrthologPseudomonas(protein_id=row['Locus Tag (Strain 1)'], strain_protein='PAO1',
                                                        ortholog_id=row['Locus Tag (Strain 2)'], strain_ortholog='PA14')
                        ortholog2 = OrthologPseudomonas(protein_id=row['Locus Tag (Strain 2)'], strain_protein='PA14',
                                                        ortholog_id=row['Locus Tag (Strain 1)'], strain_ortholog='PAO1')
                        session.add(ortholog1), session.add(ortholog2)
            session.commit()
        ortho = []

        for interaction in session.query(Interaction).all():
            interactors, ortholog_interactors = [], [[], []]
            num = 0
            for interactor in interaction.interactors:
                if interactor.type == 'protein':
                    for ortholog in interactor.pseudomonas_orthologs:
                        if ortholog is not None:
                            ortholog_interactor = session.query(Interactor).filter(Interactor.id == ortholog.ortholog_id).one()
                            ortholog_interactors[num].append(ortholog_interactor)
                else:
                    ortholog_interactors[num].append(interactor)
                num+=1

            for interactor1 in ortholog_interactors[0]:
                for interactor2 in ortholog_interactors[1]:
                    if (interactor1.type != 'metabolite') | (interactor2.type != 'metabolite'):
                        interactors.append([interactor1, interactor2])
            ortho.append(interactors)
            for interactor_pair in interactors:

                homogenous = (interactor_pair[0] == interactor_pair[1])
                ortholog_interaction = session.query(Interaction).filter(
                    (Interaction.interactors.contains(interactor_pair[0])),
                    (Interaction.interactors.contains(interactor_pair[1])),
                    (Interaction.homogenous == homogenous)).first()

                if (ortholog_interaction != None):
                    if (ortholog_interaction.ortholog_derived == None):
                        ortholog_interaction.ortholog_derived = 'confirmed from ' + interaction.strain
                    elif ('from ' + interaction.strain) not in ortholog_interaction.ortholog_derived:
                        ortholog_interaction.ortholog_derived += ', confirmed from ' + interaction.strain
                    session.commit()
                else:
                    strain = None
                    if interactor_pair[0].type == 'protein':
                        strain = interactor_pair[0].strain
                    else:
                        strain = interactor_pair[1].strain
                    ortholog_interaction = Interaction(strain=strain, type=interaction.type,
                                                       ortholog_derived='from ' + interaction.strain,
                                                       interactors=interactor_pair, homogenous = homogenous,
                                                       is_experimental = 0)
                    session.add(ortholog_interaction), session.commit()
        print(session.query(Interaction).count())

    parse_pseudomonas_orthologs()

    def parse_Ecoli():
        def parse_Ecoli_IntAct(file):
            with open(file) as csvfile:
                reader = csv.DictReader(csvfile, delimiter='\t')
                string = []
                for row in reader:
                    if (row['#ID(s) interactor A'] == '-') | (row['ID(s) interactor B'] == '-'): continue

                    string.append(row['#ID(s) interactor A'] + '$' + row['ID(s) interactor B'])

                    interactors = []
                    orthologs_A = session.query(OrthologEcoli).filter(
                        OrthologEcoli.ortholog_uniprot == row['#ID(s) interactor A'].split(':')[1]).all()
                    orthologs_B = session.query(OrthologEcoli).filter(
                        OrthologEcoli.ortholog_uniprot == row['ID(s) interactor B'].split(':')[1]).all()
                    for ortholog_A in orthologs_A:
                        if ortholog_A != None:
                            interactor_A = session.query(Interactor).filter(
                                Interactor.id == ortholog_A.protein_id).one()
                            for ortholog_B in orthologs_B:
                                if ortholog_B != None:
                                    interactor_B = session.query(Interactor).filter(
                                        Interactor.id == ortholog_B.protein_id).one()
                                    if (interactor_A.strain == interactor_B.strain):
                                        interactors.append([interactor_A, interactor_B])

                    for interactor_pair in interactors:
                        homogenous = (interactor_pair[0] == interactor_pair[1])
                        interaction = session.query(Interaction).filter(
                            (Interaction.interactors.contains(interactor_pair[0])),
                            (Interaction.interactors.contains(interactor_pair[1])),
                            (Interaction.homogenous == homogenous)).first()
                        if (interaction != None):
                            if (interaction.ortholog_derived == None):
                                interaction.ortholog_derived = 'confirmed from E.coli'
                            elif ('from E. coli' not in interaction.ortholog_derived):
                                interaction.ortholog_derived += ', confirmed from E. coli'
                            session.commit()
                        if (interaction == None):
                            interaction = Interaction(strain=interactor_pair[0].strain, interactors=interactor_pair,
                                                      type=(interactor_pair[0].type + '-' + interactor_pair[1].type),
                                                      is_experimental=0, ortholog_derived='from E. coli')
                            session.add(interaction), session.commit()
                        date = None
                        if len(row['Publication 1st author(s)'].split('(')) == 2:
                            date = row['Publication 1st author(s)'].split('(')[1][:-1]
                        reference = InteractionReference(interaction_type=row['Interaction type(s)'].split('(')[1][:-1],
                                                         interaction_id=interaction.id,
                                                         confidence_score=row['Confidence value(s)'],
                                                         pmid=row['Publication Identifier(s)'].split('med:')[1][:8],
                                                         detection_method=
                                                         row['Interaction detection method(s)'].split('(')[1][:-1],
                                                         author_last_name=row['Publication 1st author(s)'].split(' ')[
                                                             0],
                                                         publication_date=date,
                                                         publication_ref=row['Publication 1st author(s)'],
                                                         source_db=row['Source database(s)'].split('(')[1][:-1],
                                                         interactor_a_id=row['#ID(s) interactor A'].split(':')[1],
                                                         interactor_b_id=row['ID(s) interactor B'].split(':')[1])
                        session.add(reference)
                session.commit()
                print(session.query(Interaction).count())

        #parse_Ecoli_IntAct('PSICQUIC_Ecoli_IntAct.txt')

        #parse_Ecoli_IntAct('Ecoli_IntAct.txt')

        # with open('PSICQUIC_Ecoli_UniProt.txt') as csvfile:
        #     reader = csv.DictReader(csvfile, delimiter='\t')
        #     for row in reader:
        #         interactors = []
        #         orthologs_A = session.query(OrthologEcoli).filter(
        #             OrthologEcoli.ortholog_uniprot == row['#ID(s) interactor A'].split(':')[1]).all()
        #         orthologs_B = session.query(OrthologEcoli).filter(
        #             OrthologEcoli.ortholog_uniprot == row['ID(s) interactor B'].split(':')[1]).all()
        #         for ortholog_A in orthologs_A:
        #             if ortholog_A != None:
        #                 interactor_A = session.query(Interactor).filter(Interactor.id == ortholog_A.protein_id).one()
        #                 for ortholog_B in orthologs_B:
        #                     if ortholog_B != None:
        #                         interactor_B = session.query(Interactor).filter(
        #                             Interactor.id == ortholog_B.protein_id).one()
        #                         if (interactor_A.strain == interactor_B.strain):
        #                             interactors.append([interactor_A, interactor_B])
        #
        #         for interactor_pair in interactors:
        #             homogenous = (interactor_pair[0] == interactor_pair[1])
        #             interaction = session.query(Interaction).filter(
        #                 (Interaction.interactors.contains(interactor_pair[0])),
        #                 (Interaction.interactors.contains(interactor_pair[1])),
        #                 (Interaction.homogenous == homogenous)).first()
        #             if (interaction != None):
        #                 if (interaction.ortholog_derived == None):
        #                     interaction.ortholog_derived = 'confirmed from E.coli'
        #                 elif ('from E. coli' not in interaction.ortholog_derived):
        #                     interaction.ortholog_derived += ', confirmed from E. coli'
        #                 session.commit()
        #             if (interaction == None):
        #                 interaction = Interaction(strain=interactor_pair[0].strain, interactors = interactor_pair,
        #                                           type=(interactor_pair[0].type + '-' + interactor_pair[1].type),
        #                                           is_experimental=0, ortholog_derived='from E. coli')
        #                 session.add(interaction), session.commit()
        #             reference = InteractionReference(interaction_type=row['Interaction type(s)'].split('(')[1][:-1],
        #                                              interaction_id=interaction.id,
        #                                              pmid=row['Publication Identifier(s)'].split('med:')[1][:8],
        #                                              confidence_score = row['Confidence value(s)'],
        #                                              detection_method=row['Interaction detection method(s)'].split('(')[
        #                                                                   1][:-1],
        #                                              author_last_name=row['Publication 1st author(s)'].split(' ')[0],
        #                                              publication_date= row['Publication 1st author(s)'].split('(')[1][:-1],
        #                                              publication_ref=row['Publication 1st author(s)'],
        #                                              source_db=row['Source database(s)'].split('(')[1][:-1],
        #                                              interactor_a_id=row['#ID(s) interactor A'].split(':')[1],
        #                                              interactor_b_id=row['ID(s) interactor B'].split(':')[1])
        #             session.add(reference)
        #     session.commit()
        #     print(session.query(Interaction).count())
        #
        # with open('PSICQUIC_Ecoli_EBI-GOA-nonIntAct.txt') as csvfile:
        #     reader = csv.DictReader(csvfile, delimiter='\t')
        #     for row in reader:
        #         interactors = []
        #         orthologs_A = session.query(OrthologEcoli).filter(
        #             OrthologEcoli.ortholog_uniprot == row['#ID(s) interactor A'].split(':')[1]).all()
        #         orthologs_B = session.query(OrthologEcoli).filter(
        #             OrthologEcoli.ortholog_uniprot == row['ID(s) interactor B'].split(':')[1]).all()
        #         for ortholog_A in orthologs_A:
        #             if ortholog_A != None:
        #                 interactor_A = session.query(Interactor).filter(Interactor.id == ortholog_A.protein_id).one()
        #                 for ortholog_B in orthologs_B:
        #                     if ortholog_B != None:
        #                         interactor_B = session.query(Interactor).filter(
        #                             Interactor.id == ortholog_B.protein_id).one()
        #                         if (interactor_A.strain == interactor_B.strain):
        #                             interactors.append([interactor_A, interactor_B])
        #
        #         for interactor_pair in interactors:
        #             homogenous = (interactor_pair[0] == interactor_pair[1])
        #             interaction = session.query(Interaction).filter(
        #                 (Interaction.interactors.contains(interactor_pair[0])),
        #                 (Interaction.interactors.contains(interactor_pair[1])),
        #                 (Interaction.homogenous == homogenous)).first()
        #             if (interaction != None):
        #                 if (interaction.ortholog_derived == None):
        #                     interaction.ortholog_derived = 'confirmed from E.coli'
        #                 elif ('from E. coli' not in interaction.ortholog_derived):
        #                     interaction.ortholog_derived += ', confirmed from E. coli'
        #                 session.commit()
        #             if (interaction == None):
        #                 interaction = Interaction(strain=interactor_pair[0].strain, interactors=interactor_pair,
        #                                           type=(interactor_pair[0].type + '-' + interactor_pair[1].type),
        #                                           is_experimental=0, ortholog_derived='from E. coli')
        #                 session.add(interaction), session.commit()
        #             reference = InteractionReference(interaction_type=row['Interaction type(s)'].split('(')[1][:-1],
        #                                              interaction_id=interaction.id,
        #                                              pmid=row['Publication Identifier(s)'].split('med:')[1][:8],
        #                                              detection_method=row['Interaction detection method(s)'].split('(')[
        #                                                                   1][:-1],
        #                                              author_last_name=row['Publication 1st author(s)'].split(' ')[0],
        #                                              publication_date=row['Publication 1st author(s)'].split('(')[1][
        #                                                               :-1],
        #                                              publication_ref=row['Publication 1st author(s)'],
        #                                              source_db=row['Source database(s)'].split('(')[1][:-1],
        #                                              interactor_a_id=row['#ID(s) interactor A'].split(':')[1],
        #                                              interactor_b_id=row['ID(s) interactor B'].split(':')[1])
        #             session.add(reference)
        #     session.commit()
        #     print(session.query(Interaction).count())
        #
        # with open('PSICQUIC_Ecoli_IMEx.txt') as csvfile:
        #     reader = csv.DictReader(csvfile, delimiter='\t')
        #     for row in reader:
        #         if (row['#ID(s) interactor A'] == '-') | (row['ID(s) interactor B'] == '-'): continue
        #         interactors = []
        #         orthologs_A = session.query(OrthologEcoli).filter(
        #             OrthologEcoli.ortholog_uniprot == row['#ID(s) interactor A'].split(':')[1]).all()
        #         orthologs_B = session.query(OrthologEcoli).filter(
        #             OrthologEcoli.ortholog_uniprot == row['ID(s) interactor B'].split(':')[1]).all()
        #         for ortholog_A in orthologs_A:
        #             if ortholog_A != None:
        #                 interactor_A = session.query(Interactor).filter(Interactor.id == ortholog_A.protein_id).one()
        #                 for ortholog_B in orthologs_B:
        #                     if ortholog_B != None:
        #                         interactor_B = session.query(Interactor).filter(
        #                             Interactor.id == ortholog_B.protein_id).one()
        #                         if (interactor_A.strain == interactor_B.strain):
        #                             interactors.append([interactor_A, interactor_B])
        #
        #         for interactor_pair in interactors:
        #             homogenous = (interactor_pair[0] == interactor_pair[1])
        #             interaction = session.query(Interaction).filter(
        #                 (Interaction.interactors.contains(interactor_pair[0])),
        #                 (Interaction.interactors.contains(interactor_pair[1])),
        #                 (Interaction.homogenous == homogenous)).first()
        #             if (interaction != None):
        #                 if (interaction.ortholog_derived == None):
        #                     interaction.ortholog_derived = 'confirmed from E.coli'
        #                 elif ('from E. coli' not in interaction.ortholog_derived):
        #                     interaction.ortholog_derived += ', confirmed from E. coli'
        #                 session.commit()
        #             if (interaction == None):
        #                 interaction = Interaction(strain=interactor_pair[0].strain, interactors=interactor_pair,
        #                                           type=(interactor_pair[0].type + '-' + interactor_pair[1].type),
        #                                           is_experimental=0, ortholog_derived='from E. coli')
        #                 session.add(interaction), session.commit()
        #             reference = InteractionReference(interaction_type=row['Interaction type(s)'].split('(')[1][:-1],
        #                                              interaction_id=interaction.id,
        #                                              confidence_score=row['Confidence value(s)'],
        #                                              pmid=row['Publication Identifier(s)'].split('med:')[1][:8],
        #                                              detection_method=row['Interaction detection method(s)'].split('(')[
        #                                                                   1][:-1],
        #                                              author_last_name=row['Publication 1st author(s)'].split(',')[0],
        #                                              publication_date=row['Publication 1st author(s)'].split('(')[1][:-1],
        #                                              publication_ref=row['Publication 1st author(s)'],
        #                                              source_db=row['Source database(s)'].split('(')[1][:-1],
        #                                              interactor_a_id=row['#ID(s) interactor A'].split(':')[1],
        #                                              interactor_b_id=row['ID(s) interactor B'].split(':')[1])
        #             session.add(reference)
        #     session.commit()
        #     print(session.query(Interaction).count())
        # #

        # with open('PSICQUIC_Ecoli_iRefIndex.txt') as csvfile:
        #     reader = csv.DictReader(csvfile, delimiter='\t')
        #     for row in reader:
        #         if (row['#ID(s) interactor A'] == '-') | (row['ID(s) interactor B'] == '-'): continue
        #         interactors, orthologs_A, orthologs_B = [], [], []
        #         if row['#ID(s) interactor A'].split(':')[0] == 'uniprotkb':
        #             orthologs_A = session.query(OrthologEcoli).filter(
        #                 OrthologEcoli.ortholog_uniprot == row['#ID(s) interactor A'].split(':')[1]).all()
        #         elif row['#ID(s) interactor A'].split(':')[0] == 'refseq':
        #             orthologs_A = session.query(OrthologEcoli).filter(
        #                 OrthologEcoli.ortholog_refseq == row['#ID(s) interactor A'].split(':')[1]).all()
        #         if row['ID(s) interactor B'].split(':')[0] == 'uniprotkb':
        #             orthologs_B = session.query(OrthologEcoli).filter(
        #                 OrthologEcoli.ortholog_uniprot == row['ID(s) interactor B'].split(':')[1]).all()
        #         elif row['ID(s) interactor B'].split(':')[0] == 'refseq':
        #             orthologs_B = session.query(OrthologEcoli).filter(
        #                 OrthologEcoli.ortholog_refseq == row['ID(s) interactor B'].split(':')[1]).all()
        #
        #         for ortholog_A in orthologs_A:
        #             if ortholog_A != None:
        #                 interactor_A = session.query(Interactor).filter(Interactor.id == ortholog_A.protein_id).one()
        #                 for ortholog_B in orthologs_B:
        #                     if ortholog_B != None:
        #                         interactor_B = session.query(Interactor).filter(
        #                             Interactor.id == ortholog_B.protein_id).one()
        #                         if (interactor_A.strain == interactor_B.strain):
        #                             interactors.append([interactor_A, interactor_B])
        #
        #         for interactor_pair in interactors:
        #             homogenous = (interactor_pair[0] == interactor_pair[1])
        #             interaction = session.query(Interaction).filter(
        #                 (Interaction.interactors.contains(interactor_pair[0])),
        #                 (Interaction.interactors.contains(interactor_pair[1])),
        #                 (Interaction.homogenous == homogenous)).first()
        #             if (interaction != None):
        #                 if (interaction.ortholog_derived == None):
        #                     interaction.ortholog_derived = 'confirmed from E.coli'
        #                 elif ('from E. coli' not in interaction.ortholog_derived):
        #                     interaction.ortholog_derived += ', confirmed from E. coli'
        #                 session.commit()
        #             if (interaction == None):
        #                 interaction = Interaction(strain=interactor_pair[0].strain, interactors=interactor_pair,
        #                                           type=(interactor_pair[0].type + '-' + interactor_pair[1].type),
        #                                           is_experimental=0, ortholog_derived='from E. coli')
        #                 session.add(interaction), session.commit()
        #             author, date, ref, type, pmid, detection = None, None, None, None, None, None
        #             if (row['Publication 1st author(s)'] != '-'):
        #                 author = row['Publication 1st author(s)'].split('-')[0]
        #                 date = row['Publication 1st author(s)'].split('-')[1]
        #                 ref = row['Publication 1st author(s)']
        #             if (row['Interaction type(s)'] != '-'):
        #                 type = row['Interaction type(s)'].split('(')[1][:-1]
        #             if (row['Publication Identifier(s)'] != '-'):
        #                 pmid = row['Publication Identifier(s)'].split('med:')[1][:8]
        #             if (row['Interaction detection method(s)'] != '-'):
        #                 detection = row['Interaction detection method(s)'].split('(')[1][:-1]
        #             # there are more than one pmid sometimes
        #             reference = InteractionReference(interaction_type=type, interaction_id=interaction.id,
        #                                              confidence_score=row['Confidence value(s)'], pmid=pmid,
        #                                              detection_method=detection, author_last_name=author,
        #                                              publication_date=date, publication_ref=ref,
        #                                              source_db=row['Source database(s)'].split('(')[1][:-1],
        #                                              interactor_a_id=row['#ID(s) interactor A'].split(':')[1],
        #                                              interactor_b_id=row['ID(s) interactor B'].split(':')[1])
        #             session.add(reference)
        #     session.commit()
        #     print(session.query(Interaction).count())
        #
        # with open('PSICQUIC_Ecoli_mentha.txt') as csvfile:
        #     reader = csv.DictReader(csvfile, delimiter='\t')
        #     string = []
        #     for row in reader:
        #         if (row['#ID(s) interactor A'] == '-') | (row['ID(s) interactor B'] == '-'): continue
        #         string.append(row['#ID(s) interactor A'] + '$' + row['ID(s) interactor B'])
        #         interactors = []
        #         orthologs_A = session.query(OrthologEcoli).filter(
        #             OrthologEcoli.ortholog_uniprot == row['#ID(s) interactor A'].split(':')[1]).all()
        #         orthologs_B = session.query(OrthologEcoli).filter(
        #             OrthologEcoli.ortholog_uniprot == row['ID(s) interactor B'].split(':')[1]).all()
        #         for ortholog_A in orthologs_A:
        #             if ortholog_A != None:
        #                 interactor_A = session.query(Interactor).filter(Interactor.id == ortholog_A.protein_id).one()
        #                 for ortholog_B in orthologs_B:
        #                     if ortholog_B != None:
        #                         interactor_B = session.query(Interactor).filter(
        #                             Interactor.id == ortholog_B.protein_id).one()
        #                         if (interactor_A.strain == interactor_B.strain):
        #                             interactors.append([interactor_A, interactor_B])
        #
        #         for interactor_pair in interactors:
        #             homogenous = (interactor_pair[0] == interactor_pair[1])
        #             interaction = session.query(Interaction).filter(
        #                 (Interaction.interactors.contains(interactor_pair[0])),
        #                 (Interaction.interactors.contains(interactor_pair[1])),
        #                 (Interaction.homogenous == homogenous)).first()
        #             if (interaction != None):
        #                 if (interaction.ortholog_derived == None):
        #                     interaction.ortholog_derived = 'confirmed from E.coli'
        #                 elif ('from E. coli' not in interaction.ortholog_derived):
        #                     interaction.ortholog_derived += ', confirmed from E. coli'
        #                 session.commit()
        #             if (interaction == None):
        #                 interaction = Interaction(strain=interactor_pair[0].strain, interactors=interactor_pair,
        #                                           type=(interactor_pair[0].type + '-' + interactor_pair[1].type),
        #                                           is_experimental=0, ortholog_derived='from E. coli')
        #                 session.add(interaction), session.commit()
        #             reference = InteractionReference(interaction_type=row['Interaction type(s)'].split('(')[1][:-1],
        #                                              interaction_id=interaction.id,
        #                                              confidence_score=row['Confidence value(s)'],
        #                                              pmid=row['Publication Identifier(s)'].split('med:')[1][:8],
        #                                              detection_method=
        #                                              row['Interaction detection method(s)'].split('(')[1][:-1],
        #                                              source_db=row['Source database(s)'].split('(')[1][:-1])
        #             session.add(reference)
        #     session.commit()
        #     print(session.query(Interaction).count())
        #
        # with open('PSICQUIC_Ecoli_MINT.txt') as csvfile:
        #     reader = csv.DictReader(csvfile, delimiter='\t')
        #     for row in reader:
        #         if (row['#ID(s) interactor A'] == '-') | (row['ID(s) interactor B'] == '-'):
        #             continue
        #         interactors = []
        #         orthologs_A = session.query(OrthologEcoli).filter(
        #             OrthologEcoli.ortholog_uniprot == row['#ID(s) interactor A'].split(':')[1]).all()
        #         orthologs_B = session.query(OrthologEcoli).filter(
        #             OrthologEcoli.ortholog_uniprot == row['ID(s) interactor B'].split(':')[1]).all()
        #         for ortholog_A in orthologs_A:
        #             if ortholog_A != None:
        #                 interactor_A = session.query(Interactor).filter(Interactor.id == ortholog_A.protein_id).one()
        #                 for ortholog_B in orthologs_B:
        #                     if ortholog_B != None:
        #                         interactor_B = session.query(Interactor).filter(
        #                             Interactor.id == ortholog_B.protein_id).one()
        #                         if (interactor_A.strain == interactor_B.strain):
        #                             interactors.append([interactor_A, interactor_B])
        #
        #         for interactor_pair in interactors:
        #             homogenous = (interactor_pair[0] == interactor_pair[1])
        #             interaction = session.query(Interaction).filter(
        #                 (Interaction.interactors.contains(interactor_pair[0])),
        #                 (Interaction.interactors.contains(interactor_pair[1])),
        #                 (Interaction.homogenous == homogenous)).first()
        #             if (interaction != None):
        #                 if (interaction.ortholog_derived == None):
        #                     interaction.ortholog_derived = 'confirmed from E.coli'
        #                 elif ('from E. coli' not in interaction.ortholog_derived):
        #                     interaction.ortholog_derived += ', confirmed from E. coli'
        #                 session.commit()
        #             if (interaction == None):
        #                 interaction = Interaction(strain=interactor_pair[0].strain, interactors=interactor_pair,
        #                                           type=(interactor_pair[0].type + '-' + interactor_pair[1].type),
        #                                           is_experimental=0, ortholog_derived='from E. coli')
        #                 session.add(interaction), session.commit()
        #             reference = InteractionReference(interaction_type=row['Interaction type(s)'].split('(')[1][:-1],
        #                                              interaction_id=interaction.id,
        #                                              confidence_score=row['Confidence value(s)'],
        #                                              pmid=row['Publication Identifier(s)'].split('med:')[1][:8],
        #                                              detection_method=row['Interaction detection method(s)'].split('(')[
        #                                                                   1][:-1],
        #                                              author_last_name=row['Publication 1st author(s)'].split(',')[0],
        #                                              publication_date=row['Publication 1st author(s)'].split('(')[1][
        #                                                               :-1],
        #                                              publication_ref=row['Publication 1st author(s)'],
        #                                              source_db=row['Source database(s)'].split('(')[1][:-1],
        #                                              interactor_a_id=row['#ID(s) interactor A'].split(':')[1],
        #                                              interactor_b_id=row['ID(s) interactor B'].split(':')[1])
        #             session.add(reference)
        #     session.commit()
        #     print(session.query(Interaction).count())

        # with open('PSICQUIC_Ecoli_MPIDB.txt') as csvfile:
        #     reader = csv.DictReader(csvfile, delimiter='\t')
        #     for row in reader:
        #         if (row['#ID(s) interactor A'] == '-') | (row['ID(s) interactor B'] == '-'): continue
        #         interactors = []
        #         orthologs_A = session.query(OrthologEcoli).filter(
        #             OrthologEcoli.ortholog_uniprot == row['#ID(s) interactor A'].split(':')[1]).all()
        #         orthologs_B = session.query(OrthologEcoli).filter(
        #             OrthologEcoli.ortholog_uniprot == row['ID(s) interactor B'].split(':')[1]).all()
        #         for ortholog_A in orthologs_A:
        #             if ortholog_A != None:
        #                 interactor_A = session.query(Interactor).filter(Interactor.id == ortholog_A.protein_id).one()
        #                 for ortholog_B in orthologs_B:
        #                     if ortholog_B != None:
        #                         interactor_B = session.query(Interactor).filter(
        #                             Interactor.id == ortholog_B.protein_id).one()
        #                         if (interactor_A.strain == interactor_B.strain):
        #                             interactors.append([interactor_A, interactor_B])
        #
        #         for interactor_pair in interactors:
        #             homogenous = (interactor_pair[0] == interactor_pair[1])
        #             interaction = session.query(Interaction).filter(
        #                 (Interaction.interactors.contains(interactor_pair[0])),
        #                 (Interaction.interactors.contains(interactor_pair[1])),
        #                 (Interaction.homogenous == homogenous)).first()
        #             if (interaction != None):
        #                 if (interaction.ortholog_derived == None):
        #                     interaction.ortholog_derived = 'confirmed from E.coli'
        #                 elif ('from E. coli' not in interaction.ortholog_derived):
        #                     interaction.ortholog_derived += ', confirmed from E. coli'
        #                 session.commit()
        #             if (interaction == None):
        #                 interaction = Interaction(strain=interactor_pair[0].strain, interactors=interactor_pair,
        #                                           type=(interactor_pair[0].type + '-' + interactor_pair[1].type),
        #                                           is_experimental=0, ortholog_derived='from E. coli')
        #                 session.add(interaction), session.commit()
        #             reference = InteractionReference(interaction_type=row['Interaction type(s)'].split('(')[1][:-1],
        #                                              interaction_id=interaction.id,
        #                                              confidence_score=row['Confidence value(s)'],
        #                                              pmid=row['Publication Identifier(s)'].split('med:')[1][:8],
        #                                              detection_method=row['Interaction detection method(s)'].split('(')[
        #                                                                   1][:-1],
        #                                              author_last_name=row['Publication 1st author(s)'].split(' ')[0],
        #                                              publication_date=row['Publication 1st author(s)'].split('(')[1][
        #                                                               :-1],
        #                                              publication_ref=row['Publication 1st author(s)'],
        #                                              source_db=row['Source database(s)'].split('(')[1][:-1],
        #                                              interactor_a_id=row['#ID(s) interactor A'].split(':')[1],
        #                                              interactor_b_id=row['ID(s) interactor B'].split(':')[1])
        #             session.add(reference)
        #     session.commit()
        #     print(session.query(Interaction).count())
        #
        # with open('Ecoli_DIP.txt') as csvfile:
        #     reader = csv.DictReader(csvfile, delimiter='\t')
        #     for row in reader:
        #         interactors = []
        #
        #         ids_A = row['ID interactor A'].split('|')
        #         ids_B = row['ID interactor B'].split('|')
        #         refseq_A, uniprotkb_A, refseq_B, uniprotkb_B = '', '', '', ''
        #         for id in ids_A:
        #             fields = id.split(':')
        #             if (fields[0] == 'refseq'):
        #                 refseq_A = fields[1]
        #             elif (fields[0] == 'uniprotkb'):
        #                 uniprotkb_A = fields[1]
        #         for id in ids_B:
        #             fields = id.split(':')
        #             if (fields[0] == 'refseq'):
        #                 refseq_B = fields[1]
        #             elif (fields[0] == 'uniprotkb'):
        #                 uniprotkb_B = fields[1]
        #
        #         orthologs_A, orthologs_B = [], []
        #         if (uniprotkb_A != ''):
        #             orthologs_A = session.query(OrthologEcoli).filter(OrthologEcoli.ortholog_uniprot == uniprotkb_A).all()
        #         if ((len(orthologs_A) == 0) & (refseq_A != '')):
        #             orthologs_A = session.query(OrthologEcoli).filter(OrthologEcoli.ortholog_refseq == refseq_A).all()
        #         if (uniprotkb_B != ''):
        #             orthologs_B = session.query(OrthologEcoli).filter(OrthologEcoli.ortholog_uniprot == uniprotkb_B).all()
        #         if ((len(orthologs_B) == 0) & (refseq_B != '')):
        #             orthologs_B = session.query(OrthologEcoli).filter(OrthologEcoli.ortholog_refseq == refseq_B).all()
        #
        #         for ortholog_A in orthologs_A:
        #             if ortholog_A != None:
        #                 interactor_A = session.query(Interactor).filter(Interactor.id == ortholog_A.protein_id).one()
        #                 for ortholog_B in orthologs_B:
        #                     if ortholog_B != None:
        #                         interactor_B = session.query(Interactor).filter(
        #                             Interactor.id == ortholog_B.protein_id).one()
        #                         if (interactor_A.strain == interactor_B.strain):
        #                             interactors.append([interactor_A, interactor_B])
        #
        #         for interactor_pair in interactors:
        #             homogenous = (interactor_pair[0] == interactor_pair[1])
        #             interaction = session.query(Interaction).filter(
        #                 (Interaction.interactors.contains(interactor_pair[0])),
        #                 (Interaction.interactors.contains(interactor_pair[1])),
        #                 (Interaction.homogenous == homogenous)).first()
        #             if (interaction != None):
        #                 if (interaction.ortholog_derived == None):
        #                     interaction.ortholog_derived = 'confirmed from E.coli'
        #                 elif ('from E. coli' not in interaction.ortholog_derived):
        #                     interaction.ortholog_derived += ', confirmed from E. coli'
        #                 session.commit()
        #             if (interaction == None):
        #                 interaction = Interaction(strain=interactor_pair[0].strain, interactors=interactor_pair,
        #                                           type=(interactor_pair[0].type + '-' + interactor_pair[1].type),
        #                                           is_experimental=0, ortholog_derived='from E. coli')
        #                 session.add(interaction), session.commit()
        #             detections, pmids, types, list = [], [], [], []
        #             if (row['Interaction detection method(s)'] != '-'):
        #                 detections = row['Interaction detection method(s)'].split('|')
        #                 list.append(detections)
        #             if (row['Publication Identifier(s)'] != '-'):
        #                 pmids = row['Publication Identifier(s)'].split('|')
        #                 list.append(pmids)
        #             if (row['Interaction type(s)']):
        #                 types = row['Interaction type(s)'].split('|')
        #                 list.append(types)
        #
        #             if (len(list) != 0) & (all((len(item) == len(list[0])) for item in list)):
        #                 for num in range(0, len(list[0])):
        #                     type, pmid, detection = None, None, None
        #                     for item in list:
        #                         if (item == types):
        #                             type = types[num].split('(')[1][:-1]
        #                         if (item == pmids):
        #                             pmid = pmids[num*2].split('med:')[1][:8]
        #                         if (item == detections):
        #                             detection = detections[num].split('(')[1][:-1]
        #                     # there are more than one pmid sometimes
        #                     reference = InteractionReference(interaction_type=type, interaction_id=interaction.id,
        #                                                      pmid=pmid, detection_method=detection,
        #                                                      source_db=row['Source database(s)'].split('(')[1][:-1])
        #                     session.add(reference)
        #     session.commit()
        #     print(session.query(Interaction).count())

        # with open('Ecoli_RegulonDB.csv') as csvfile:
        #     reader = csv.DictReader(csvfile)
        #     for row in reader:
        #         interactors = []
        #
        #         orthologs_A = session.query(OrthologEcoli).filter(
        #             OrthologEcoli.ortholog_name == (row['TF name'][0].lower() + row['TF name'][1:])).all()
        #         orthologs_B = session.query(OrthologEcoli).filter(
        #             OrthologEcoli.ortholog_name == row['Regulated gene']).all()
        #
        #         for ortholog_A in orthologs_A:
        #             if ortholog_A != None:
        #                 interactor_A = session.query(Interactor).filter(Interactor.id == ortholog_A.protein_id).one()
        #                 for ortholog_B in orthologs_B:
        #                     if ortholog_B != None:
        #                         interactor_B = session.query(Interactor).filter(
        #                             Interactor.id == ortholog_B.protein_id).one()
        #                         if (interactor_A.strain == interactor_B.strain):
        #                             interactors.append([interactor_A, interactor_B])
        #
        #         for interactor_pair in interactors:
        #             homogenous = (interactor_pair[0] == interactor_pair[1])
        #             interaction = session.query(Interaction).filter(
        #                 (Interaction.interactors.contains(interactor_pair[0])),
        #                 (Interaction.interactors.contains(interactor_pair[1])),
        #                 (Interaction.homogenous == homogenous)).first()
        #             if (interaction != None):
        #                 if (interaction.ortholog_derived == None):
        #                     interaction.ortholog_derived = 'confirmed from E.coli'
        #                 elif ('from E. coli' not in interaction.ortholog_derived):
        #                     interaction.ortholog_derived += ', confirmed from E. coli'
        #                 session.commit()
        #             if (interaction == None):
        #                 interaction = Interaction(strain=interactor_pair[0].strain, interactors=interactor_pair,
        #                                           type=(interactor_pair[0].type + '-' + interactor_pair[1].type),
        #                                           is_experimental=0, ortholog_derived='from E. coli')
        #                 session.add(interaction), session.commit()
        #             reference = InteractionReference(interaction_id=interaction.id,
        #                                              interaction_type='TF/sigma-binding site (' +
        #                                                               row['Regulatory effect'] + 'regulation)',
        #                                              detection_method=row['Evidence'],
        #                                              interaction_full_name=row['TF name'].lower() + ' regulates(' +
        #                                                                    row['Regulatory effect'] + ') ' +
        #                                                                    row['Regulated gene'],
        #                                              source_db = 'RegulonDB')
        #             session.add(reference)
        #     session.commit()
        #     print(session.query(InteractionReference).count())


    def get_ecocyc_compounds():
        string = ''
        with open("ECOCYC_paths.txt") as csvfile:
            path_reader = csv.DictReader(csvfile)
            for path_row in path_reader:
                interactor_file_name = "ecocyc_sif_files/" + path_row["Pathways"]  + "_interactors.txt"
                if not exists(interactor_file_name): continue
                with open(interactor_file_name) as pathfile:
                    interactor_reader = csv.DictReader(pathfile)
                    for interactor_row in interactor_reader:
                        if interactor_row["PARTICIPANT_TYPE"] != "SmallMoleculeReference": continue
                        xrefs = interactor_row["UNIFICATION_XREF"]
                        name = interactor_row['PARTICIPANT']
                        if (name not in ecocyc_compounds):
                            ecocyc_compounds[name] = {}
                            id_list = [['PubChem-compound:', 'PubChem'], ['ChEBI:', 'ChEBI'], ['CAS:', 'CAS'],
                                       ['KEGG LIGAND:', 'KEGG'], ['EcoCyc:', 'EcoCyc']]
                            for id_type in id_list:
                                ecocyc_compounds[name][id_type[1]] = None
                                if (len(xrefs.split(id_type[0])) < 2): continue
                                ecocyc_compounds[name][id_type[1]] = xrefs.split(id_type[0])[1].split(';')[0]
        print(string)
        # for metabolite in ecocyc_compounds:
        #     if 'KEGG' in metabolite:
        #         if (session.query(Metabolite).filter(Metabolite.id == metabolite["KEGG"]).first() is not None):
        #             interactor = session.query(Metabolite.filter(Metabolite.id == metabolite["KEGG"])).one()
        #             interactor.CAS = metabolite['CAS']
        #             interactor.ChEBI = metabolite['ChEBI']
        #             interactor.EcoCyc = metabolite['EcoCyc']
        #             interactor.id = metabolite['EcoCyc']


    def parse_ecocyc(strain):
        get_ecocyc_compounds()
        with open("ECOCYC_paths.txt") as csvfile:
            path_reader = csv.DictReader(csvfile)
            for path_row in path_reader:
                interaction_file_name = "ecocyc_sif_files/" + path_row["Pathways"] + "_interactions.txt"
                if not exists(interaction_file_name): continue
                with open(interaction_file_name) as pathfile:
                    interaction_reader = csv.DictReader(pathfile)
                    for interaction_row in interaction_reader:
                        interactors = [[],[]]
                        new_metabolites = [[],[]]

                        A_id = interaction_row['PARTICIPANT_A']
                        B_id = interaction_row['PARTICIPANT_B']

                        if (A_id not in ecocyc_compounds):
                            for ortholog in (session.query(OrthologEcoli).filter(
                                    (OrthologEcoli.ortholog_uniprot == A_id),
                                    (OrthologEcoli.strain_protein == strain)).all()):
                                if (ortholog is not None): interactors[0].append(ortholog.protein)
                        else:
                            A_ecocyc = ecocyc_compounds[A_id]["EcoCyc"]
                            if (session.query(Metabolite).filter(Metabolite.id == A_ecocyc).first() is not None):
                                interactors[0].append(session.query(Metabolite).filter(
                                    Metabolite.id == A_ecocyc).one())
                            else:
                                new_metabolites[0].append(A_ecocyc)

                        if (B_id not in ecocyc_compounds):
                            for ortholog in (session.query(OrthologEcoli).filter(
                                    (OrthologEcoli.ortholog_uniprot == B_id),
                                    (OrthologEcoli.strain_protein == strain)).all()):
                                if (ortholog is not None): interactors[1].append(ortholog.protein)
                        else:
                            B_ecocyc = ecocyc_compounds[B_id]["EcoCyc"]
                            if (session.query(Metabolite).filter(Metabolite.id == B_ecocyc).first() is not None):
                                interactors[1].append(session.query(Metabolite).filter(
                                    Metabolite.id == B_ecocyc).one())
                            else:
                                new_metabolites[1].append(B_ecocyc)
                                print(B_id, B_ecocyc, path_row['Pathways'])

                        interactor_pairs = []
                        for interactor1 in interactors[0]:
                            for interactor2 in interactors[1]:
                                if (interactor1.type != 'metabolite') | (interactor2.type != 'metabolite'):
                                    interactor_pairs.append([interactor1, interactor2])

                        for interactor1 in interactors[0]:
                            for interactor2 in new_metabolites[1]:
                                if (interactor1.type != 'metabolite'):
                                    metabolite = session.query(Metabolite).filter(Metabolite.id == interactor2).first()
                                    if (metabolite is None):
                                        metabolite = Metabolite(id = interactor2, EcoCyc = interactor2,
                                                                PubChem = ecocyc_compounds[B_id]['PubChem'],
                                                                KEGG = ecocyc_compounds[B_id]['KEGG'],
                                                                CAS = ecocyc_compounds[B_id]['CAS'],
                                                                ChEBI = ecocyc_compounds[B_id]['ChEBI'], name = B_id)
                                        session.add(metabolite), session.commit()
                                    interactor_pairs.append([interactor1, metabolite])

                        for interactor1 in interactors[1]:
                            for interactor2 in new_metabolites[0]:
                                if (interactor1.type != 'metabolite'):
                                    metabolite = session.query(Metabolite).filter(Metabolite.id == interactor2).first()
                                    if (metabolite is None):
                                        metabolite = Metabolite(id=interactor2, EcoCyc=interactor2,
                                                                PubChem=ecocyc_compounds[A_id]['PubChem'],
                                                                KEGG=ecocyc_compounds[A_id]['KEGG'],
                                                                CAS=ecocyc_compounds[A_id]['CAS'],
                                                                ChEBI = ecocyc_compounds[A_id]['ChEBI'], name = A_id)
                                        session.add(metabolite), session.commit()
                                    interactor_pairs.append([metabolite, interactor1])

                        if (len(interactor_pairs) == 0): continue

                        for interactor_pair in interactor_pairs:
                            homogenous = (interactor_pair[0] == interactor_pair[1])
                            interaction = session.query(Interaction).filter(
                                (Interaction.interactors.contains(interactor_pair[0])),
                                (Interaction.interactors.contains(interactor_pair[1])),
                                (Interaction.homogenous == homogenous)).first()

                            if (interaction == None):
                                interaction = Interaction(type=interactor_pair[0].type + '-' + interactor_pair[1].type,
                                                          strain=strain, homogenous=homogenous,
                                                          interactors=interactor_pair,
                                                          comment= interactor_pair[0].id +
                                                                   interaction_row["INTERACTION_TYPE"] +
                                                                   interactor_pair[1].id,
                                                          ortholog_derived = 'from E. coli')
                                session.add(interaction), session.commit()
                            else:
                                interaction = session.query(Interaction).filter(
                                    (Interaction.interactors.contains(interactor_pair[0])),
                                    (Interaction.interactors.contains(interactor_pair[1])),
                                    (Interaction.homogenous == homogenous)).one()
                                if (interaction.ortholog_derived == None):
                                    interaction.ortholog_derived = 'confirmed from E. coli'
                                else:
                                    if ('from E.coli' not in interaction.ortholog_derived):
                                        interaction.ortholog_derived += ', confirmed from E.coli'

                            for pmid in interaction_row["INTERACTION_PUBMED_ID"].split(';'):
                                if (session.query(InteractionReference).filter(
                                        InteractionReference.interaction_id == interaction.id,
                                        InteractionReference.pmid == pmid).first() is None):
                                    reference = InteractionReference(interaction_id = interaction.id, pmid = pmid)
                                    session.add(reference), session.commit()

                            if (session.query(InteractionXref).filter(
                                    (InteractionXref.interaction_id == interaction.id),
                                    (InteractionXref.data_source == 'EcoCyc')).first() == None):
                                xref = InteractionXref(interaction_id=interaction.id, data_source='EcoCyc')
                                session.add(xref), session.commit()
                            session.commit()
        print(session.query(Interaction).count())


    parse_kegg('eco', 'PAO1')
    parse_kegg('eco', 'PA14')
    parse_ecocyc('PAO1')
    parse_ecocyc('PA14')
    parse_Ecoli()
    print([session.query(Interaction).filter(Interaction.strain == 'PAO1').count(),
           session.query(Interaction).filter(Interaction.strain == 'PA14').count()])