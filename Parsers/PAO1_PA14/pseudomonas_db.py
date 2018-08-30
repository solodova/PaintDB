import csv
from Schema1 import Protein, ProteinComplex, Interactor, InteractorXref, Localization, GeneOntology

def parse_pseudomonasdb(session):
    get_gene_info('Data/PAO1/Gene_Info/pseudomonas_db_info.csv', 'PAO1', session)
    #get_xrefs('Data/PAO1/Gene_Info/pseudomonas_db_xrefs.csv', 'PAO1', session)
    #get_localizations('Data/PAO1/Gene_Info/pseudomonas_db_localizations.csv', session)
    #get_ontology('Data/PAO1/Gene_Info/pseudomonas_db_go.csv', session)

    get_gene_info('Data/PA14/Gene_Info/pseudomonas_db_info.csv', 'PA14', session)
    #get_xrefs('Data/PA14/Gene_Info/pseudomonas_db_xrefs.csv', 'PA14', session)
    # get_localizations('Data/PA14/Gene_Info/pseudomonas_db_localizations.csv', session)
    # get_ontology('Data/PA14/Gene_Info/pseudomonas_db_go.csv', session)


def get_gene_info(file, strain, session):
    with open(file) as csvfile:
        # ignore all other fields
        # need to set field names because first row is not the column names
        fieldnames = ["Sequence", "Locus Tag", "Feature Type", "Start", "End", "Strand", "Name", "Product Name",
                      "Synonyms", "NCBI Accession"]
        reader = csv.DictReader(csvfile, fieldnames=fieldnames)
        row_num = 0
        for row in reader:
            # skip first three rows since they don't contain interactor info
            if (row['Feature Type'] == 'CDS') and (row_num >= 3):
                # trim trailing " character from locus tag
                protein = Protein(id=row['Locus Tag'][:-1], name=row['Name'], type='p', strain=strain,
                                  product_name=row['Product Name'], ncbi_acc=row['NCBI Accession'])
                session.add(protein)
            row_num += 1
        session.commit()


def get_xrefs(file, strain, session):
    with open(file) as csvfile:
        # types = ['RefSeq Accession', 'UniProtKB Accession', 'UniProtKB ID', 'GI Number', 'Uniparc', 'UniRef100 ID',
        #          'UniRef90 ID', 'UniRef50 ID']
        types = ['UniProtKB Accession', 'RefSeq Accession']
        reader = csv.DictReader(csvfile)
        for row in reader:
            if (row['Feature Type'] == 'CDS') and \
                    (session.query(Interactor).filter(Interactor.id == row['Locus Tag']).first() is not None):
                for type in types:
                    # if no accession present, don't include xref
                    if row[type] == '': continue
                    # create and add xref if it isn't present already (need to check if xref already exists since
                    # there are duplicate xrefs in file
                    if session.query(InteractorXref).filter(InteractorXref.interactor_id == row['Locus Tag'],
                                                             InteractorXref.accession == row[type]).first() is None:
                        xref = InteractorXref(accession=row[type], interactor_id=row['Locus Tag'], source=type)
                        session.add(xref)

                    # for uniprotkb ids, also need to set uniprotkb attribute of corresponding protein
                    if type == 'UniProtKB Accession':
                        # for uniprotkb ids, also need to set uniprotkb attribute of corresponding protein
                        interactor = session.query(Interactor).filter(Interactor.id == row['Locus Tag']).one()
                        interactor.uniprotkb = row[type]
                        session.commit()

                        num_interactors = session.query(Protein).filter(Protein.uniprotkb == row[type]).count()
                        # if two interactors have the same uniprotkb id, they are part of a protein complex together
                        # create a new interactor with id being the uniprotkb id, and 'protein complex' as product name
                        if num_interactors == 2:
                            interactor = ProteinComplex(id=row[type], strain=strain, type='pc')
                            session.add(interactor)
                session.commit()


def get_localizations(file, session):
    with open(file) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if session.query(Interactor).filter(Interactor.id == row['Locus Tag']).first() is not None:
                localization = Localization(protein_id=row['Locus Tag'], localization=row['Subcellular Localization'],
                                            confidence=row['Confidence'])
                session.add(localization)
        session.commit()


def get_ontology(file, session):
    with open(file) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if session.query(Interactor).filter(Interactor.id == row['Locus Tag']).first() is not None:
                ontology = GeneOntology(gene_id=row['Locus Tag'], accession=row['Accession'], go_term=row['GO Term'],
                                        evidence_code=row['GO Evidence Code'], pmid=row['PMID'],
                                        eco_code=row['Evidence Ontology ECO Code'],
                                        eco_term=row['Evidence Ontology Term'])
                session.add(ontology)
        session.commit()
