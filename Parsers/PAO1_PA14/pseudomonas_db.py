import csv
from DB_schema import Protein, Interactor, Localization, GeneOntology

def parse(session):
    get_gene_info('Data/PAO1/Gene_Info/pseudomonas_db_info.csv', 'PAO1', session)
    get_uniprot('Data/PAO1/Gene_Info/pseudomonas_db_xrefs.csv', 'PAO1', session)
    get_localizations('Data/PAO1/Gene_Info/pseudomonas_db_localizations.csv', session)
    get_ontology('Data/PAO1/Gene_Info/pseudomonas_db_go.csv', session)

    get_gene_info('Data/PA14/Gene_Info/pseudomonas_db_info.csv', 'PA14', session)
    get_uniprot('Data/PA14/Gene_Info/pseudomonas_db_xrefs.csv', 'PA14', session)
    get_localizations('Data/PA14/Gene_Info/pseudomonas_db_localizations.csv', session)
    get_ontology('Data/PA14/Gene_Info/pseudomonas_db_go.csv', session)


def get_gene_info(file, strain, session):
    with open(file) as csvfile:
        # ignore all other fields
        # need to set field names because first row is not the column names
        fieldnames = ["Sequence", "Locus Tag", "Feature Type", "Start", "End", "Strand", "Name", "Product Name",
                      "Synonyms", "NCBI Accession"]
        reader = csv.DictReader(csvfile, fieldnames=fieldnames)
        row_num = 0
        proteins = []
        for row in reader:
            # skip first three rows since they don't contain interactor info
            if (row['Feature Type'] == 'CDS') and (row_num >= 3):
                ncbi_acc, name = None, None
                if row['NCBI Accession'] != '':
                    ncbi_acc = row['NCBI Accession']
                if row['Name'] != '':
                    name = row['Name']
                # trim trailing " character from locus tag
                proteins.append(Protein(id=row['Locus Tag'][:-1], name=name, strain=strain,
                                        product_name=row['Product Name'], ncbi_acc=ncbi_acc))
            row_num += 1
        session.add_all(proteins)
    session.commit()
    print('gene info')


def get_uniprot(file, strain, session):
    with open(file) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if row['UniProtKB Accession'] == '': continue
            if (row['Feature Type'] == 'CDS') and \
                    (session.query(Interactor).get(row['Locus Tag']) is not None):
                # if no accession present, don't include xref
                interactor = session.query(Interactor).get(row['Locus Tag'])
                interactor.uniprotkb = row['UniProtKB Accession']
                num_interactors = session.query(Protein).filter_by(uniprotkb = row['UniProtKB Accession']).count()
                # if two interactors have the same uniprotkb id, they are part of a protein complex together
                # create a new interactor with id being the uniprotkb id, and 'protein complex' as product name
                if num_interactors == 2:
                    interactor = Protein(id=row['UniProtKB Accession'], strain=strain, uniprotkb = 'pc')
                    session.add(interactor)
    session.commit()
    print('xrefs')


# function to collect all protein localizations from file
def get_localizations(file, session):
    with open(file) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            interactor = session.query(Interactor).get(row['Locus Tag'])
            # only add localization if the interactor exists in the database
            if interactor is not None:
                # check if localization already exists
                localization = session.query(Localization).filter_by(localization=row['Subcellular Localization'],
                                                                     confidence=row['Confidence']).first()
                # create localization and add it to interactor localization list if it doesnt exist yet
                if localization is None:
                    interactor.localizations.append(Localization(localization=row['Subcellular Localization'],
                                                                 confidence=row['Confidence']))
                # if the localization already exists, check if it's listed in the interactors localizations,
                # if not, add it
                elif localization not in interactor.localizations:
                    interactor.localizations.append(localization)
    session.commit()
    print('localizations')

# function to get all ontologies for proteins from file
def get_ontology(file, session):
    with open(file) as csvfile:
        reader = csv.DictReader(csvfile)
        ontologies = []
        for row in reader:
            # only create ontology if interactor exists in the database
            if session.query(Interactor).get(row['Locus Tag']) is not None:
                pmid, evidence_code, eco_code = None, None, None
                if row['PMID'] != '':
                    pmid=row['PMID']
                if row['GO Evidence Code'] != '':
                    evidence_code = row['GO Evidence Code']
                if row['Evidence Ontology ECO Code'] != '':
                    eco_code = row['Evidence Ontology ECO Code']
                ontologies.append(GeneOntology(gene_id=row['Locus Tag'], accession=row['Accession'],
                                               go_term=row['GO Term'], evidence_code=evidence_code, pmid=pmid,
                                               eco_code=eco_code, eco_term=row['Evidence Ontology Term']))
        session.add_all(ontologies)
    session.commit()
    print('ontologies')
