import csv
from Schema1 import Protein, Interactor, ProteinXref, Localization, GeneOntology

def parse_pseudomonasdb(session):
    getGeneInfo('PAO1/Gene_Info/pseudomonas_db_info.csv', 'PAO1', session)
    getGeneInfo('PA14/Gene_Info/pseudomonas_db_info.csv', 'PA14', session)

    getXrefs('PAO1/Gene_Info/pseudomonas_db_xrefs.csv', 'PAO1', session)
    getXrefs('PA14/Gene_Info/pseudomonas_db_xrefs.csv', 'PA14', session)

    #getLocalizations('PAO1/Gene_Info/pseudomonas_db_localizations.csv', session)
    # getLocalizations('PA14/Gene_Info/pseudomonas_db_localizations.csv', session)

    # getOntology('PAO1/Gene_Info/pseudomonas_db_GO.csv', session)
    # getOntology('PA14/Gene_Info/pseudomonas_db_GO.csv', session)

def getGeneInfo(file, strain, session):
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


def getXrefs(file, strain, session):
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


def getLocalizations(file, session):
    with open(file) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            localization = Localization(protein_id=row['Locus Tag'], localization=row['Subcellular Localization'],
                                        confidence=row['Confidence'])
            session.add(localization)
        session.commit()


def getOntology(file, session):
    with open(file) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            ontology = GeneOntology(protein_id=row['Locus Tag'], accession=row['Accession'], go_term=row['GO Term'],
                                    evidence_code=row['GO Evidence Code'], pmid=row['PMID'],
                                    eco_code=row['Evidence Ontology ECO Code'],
                                    eco_term=row['Evidence Ontology Term'])
            session.add(ontology)
        session.commit()
