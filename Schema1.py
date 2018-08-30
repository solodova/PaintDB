import sqlalchemy
from sqlalchemy import create_engine, Column, Integer, String, ForeignKey, DateTime, Table, Float
from sqlalchemy.orm import relationship
from sqlalchemy.orm import backref
from sqlalchemy.ext.declarative import declarative_base

# stopped at chromatography tech
experimental_psimi = ['MI:0045', 'MI:0401', 'MI:0400','MI:0004', 'MI:2288', 'MI:0019', 'MI:0006', 'MI:0007',
                          'MI:0402', 'MI:0225', 'MI:1218', 'MI:1028', 'MI:1029', 'MI:0858', 'MI:0946', 'MI:1017',
                          'MI:0729', 'MI:0096', 'MI:0963', 'MI:0676', 'MI:0008', 'MI:0225', 'MI:0081', 'MI:0089',
                          'MI:0678', 'MI:0695', 'MI:0092', 'MI:0946', 'MI:0095', 'MI:0921', 'MI:0405', 'MI:0034',
                          'MI:0009', 'MI:0073', 'MI:0084', 'MI:0048', 'MI:0899', 'MI:0900', 'MI:0066', 'MI:0108',
                          'MI:0098', 'MI:0115', 'MI:1087', 'MI:1031', 'MI:0813', 'MI:0440', 'MI:0892', 'MI:2189',
                          'MI:0947', 'MI:0411', 'MI:0047', 'MI:2283', 'MI:0049', 'MI:2167', 'MI:0657', 'MI:1232']

def is_experimental_psimi(psi_code):
    return psi_code in experimental_psimi

Base = declarative_base()

protein_reference = Table('protein_reference', Base.metadata,
                          Column('protein_id', String, ForeignKey('protein.id')),
                          Column('reference_pmid', String, ForeignKey('reference.pmid')))

interaction_participant = Table('interaction_participant', Base.metadata,
                                Column('interactor_id', String, ForeignKey('interactor.id')),
                                Column('interaction_id', Integer, ForeignKey('interaction.id')))

interaction_references = Table('interaction_references', Base.metadata,
                                Column('interaction_id', String, ForeignKey('interaction.id')),
                                Column('reference_id', Integer, ForeignKey('interaction_reference.id')))

interaction_sources = Table('interaction_sources', Base.metadata,
                                Column('interaction_id', String, ForeignKey('interaction.id')),
                                Column('data_source', Integer, ForeignKey('interaction_source.data_source')))


class Interactor(Base):
    __tablename__ = 'interactor'

    # protein: locus id (eg. PA0001)
    # protein complex: uniprotkb id (eg. Q9HT76)
    # metabolite: kegg id (eg. )
    id = Column(String, primary_key=True)
    # protein: eg. "", pyrE
    # protein complex: None
    # metabolite:
    name = Column(String)
    type = Column(String)

    __mapper_args__ = {
        'polymorphic_identity': 'interactor',
        'polymorphic_on': type
    }

    interactions = relationship("Interaction", secondary=interaction_participant, backref="interactors")
    xrefs = relationship("InteractorXref", backref = 'interactor')


class Metabolite(Interactor):
    __tablename__ = 'metabolite'

    #KEGG, then EcoCyc code, then pubchem, then chebi
    id = Column(String, ForeignKey('interactor.id'), primary_key=True)
    # eg.
    kegg = Column(String)
    # eg.
    pubchem = Column(String)
    # eg.
    cas = Column(String)
    # eg.
    chebi = Column(String)
    # eg.
    ecocyc = Column(String)

    __mapper_args__ = {
        'polymorphic_identity': 'm',
    }


class Protein(Interactor):
    __tablename__ = 'protein'

    # locus id (eg. PA0001)
    id = Column(String, ForeignKey('interactor.id'), primary_key=True)
    strain = Column(String)
    # monomeric protein: product name (eg. orotate phosphoribosyltransferase, hypothetical protein)
    product_name = Column(String)
    # eg. NP_254184.1
    ncbi_acc = Column(String)
    # eg. Q9HT76
    uniprotkb = Column(String)

    __mapper_args__ = {
        'polymorphic_identity': 'p',
    }

    localizations = relationship("Localization", backref="protein")
    references = relationship("Reference", secondary=protein_reference, backref="proteins")
    pseudomonas_orthologs = relationship("OrthologPseudomonas", backref ="protein")
    ecoli_ortholgs = relationship("OrthologEcoli", backref = "protein")


class ProteinComplex(Interactor):
    __tablename__ = 'protein_complex'

    # uniprotkb id (eg. Q9HT76)
    id = Column(String, ForeignKey('interactor.id'), primary_key = True)
    strain = Column(String)

    __mapper_args__ = {
        'polymorphic_identity': 'pc',
    }

class InteractorXref(Base):
    # note that some PAO1 proteins will have multiple xrefs from the same source
    __tablename__ = 'interactor_xref'

    interactor_id = Column(String, ForeignKey("interactor.id"), primary_key=True)
    # proteins:
    #   - RefSeq Accession: NP_254184.1
    #   - UniProtKB Accession: Q9HT76
    #   - UniProtKB ID: Q9HT76_PSEAE
    #   - GI Number: 15600690
    #   - Uniparc: UPI00000C6038
    #   - UniRef100 ID: UniRef100_Q9HT76
    #   - UniRef90 ID: UniRef90_A6VEX4
    #   - UniRef50 ID: UniRef50_U2E6Q5
    accession = Column(String, primary_key=True)
    source = Column(String)


class Reference(Base):
    __tablename__ = 'reference'

    pmid = Column(String, primary_key=True)
    pub_date = Column(String)
    author_ln = Column(String)
    full_ref = Column(String)


class GeneOntology(Base):
    __tablename__ = 'gene_ontology'

    # automatically generated
    id = Column(Integer, primary_key=True, autoincrement=True)
    gene_id = Column(String, ForeignKey('protein.id'))
    # GO accession (eg. GO:0005524)
    accession = Column(String)
    # eg. ATP binding
    go_term = Column(String)
    # eg. ISM
    evidence_code = Column(String)
    # eg. ECO:0000259
    eco_code = Column(String)
    # eg. match to InterPro signature evidence used in automatic assertion
    eco_term = Column(String)
    # eg PMID:24451626
    pmid = Column(String)


class Localization(Base):
    __tablename__ = 'localization'

    id = Column(Integer, primary_key=True, autoincrement=True)
    # eg. PA0001
    protein_id = Column(String, ForeignKey('protein.id'))
    # eg. Cytoplasmic
    localization = Column(String)
    # eg. Class 3
    confidence = Column(String)


class OrthologPseudomonas(Base):
    __tablename__ = 'ortholog_pseudomonas'

    protein_id  = Column(String, ForeignKey('protein.id'), primary_key=True)
    ortholog_id = Column(String, primary_key=True)
    strain_protein = Column(String, primary_key=True)
    strain_ortholog = Column(String)


class OrthologEcoli(Base):
    __tablename__ = 'ortholog_ecoli'

    protein_id = Column(String, ForeignKey('protein.id'), primary_key=True)
    ortholog_id = Column(String, primary_key=True)
    ortholog_uniprot = Column(String)
    ortholog_refseq = Column(String)
    ortholog_name = Column(String)
    strain_protein = Column(String)
    ortholuge_classification = Column(String)


class Interaction(Base):
    __tablename__ = 'interaction'

    id = Column(Integer, primary_key=True, autoincrement=True)
    strain = Column(String)
    # remains None unless definitively determined
    #is_experimental = Column(Integer)
    # remains None unless confirmation from ortholog or derived from ortholog
    ortholog_derived = Column(String)
    type = Column(String)
    # remains None if the interaction is not tfbs, otherwise describes interaction
    is_tfbs = Column(Integer)
    # whether the two interactors are the same (0 or 1)
    homogenous = Column(Integer)

    xrefs = relationship("InteractionXref", backref="interaction")
    references = relationship("InteractionReference", secondary=interaction_references, backref="interactions")
    sources = relationship("InteractionSource", secondary=interaction_sources, backref="interactions")


class InteractionReference(Base):
    __tablename__ = 'interaction_reference'
    id = Column(Integer, primary_key=True, autoincrement=True)
    #interaction_id = Column(Integer, ForeignKey('interaction.id'))

    # PSI MI term (except Geoff)
    psimi_detection = Column(String)
    detection_method = Column(String)
    author_ln = Column(String)
    pub_date = Column(String)
    pmid = Column(String)
    psimi_type = Column(String)
    # remains None if '-', otherwise psimi OR just 'predicted'
    interaction_type = Column(String)
    psimi_db = Column(String)
    source_db = Column(String)
    # when there is additional comment about the interaction (eg. more detailed description)
    comment=Column(String)
    # includes type of confidence score, usually preceding confidence value with type/source:value
    confidence = Column(String)
    interactor_a = Column(String)
    interactor_b = Column(String)



class InteractionXref(Base):
    __tablename__ = 'interaction_xref'

    accession = Column(String, primary_key=True)
    interaction_id = Column(Integer, ForeignKey('interaction.id'), primary_key=True)
    data_source = Column(String)

class InteractionSource(Base):
    __tablename__ ='interaction_source'

    data_source=Column(String, primary_key=True)
    is_experimental = Column(Integer, primary_key=True)






