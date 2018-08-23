import sqlalchemy
from sqlalchemy import create_engine, Column, Integer, String, ForeignKey, DateTime, Table, Float
from sqlalchemy.orm import relationship
from sqlalchemy.orm import backref
from sqlalchemy.ext.declarative import declarative_base


Base = declarative_base()

protein_reference = Table('protein_reference', Base.metadata,
                             Column('protein_id', String, ForeignKey('protein.id')),
                             Column('reference_pmid', String, ForeignKey('reference.pmid')))

interaction_participant = Table('interaction_participant', Base.metadata,
                                Column('interactor_id', String, ForeignKey('interactor.id')),
                                Column('interaction_id', Integer, ForeignKey('interaction.id')))


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
    xrefs = relationship("InteractorXrefs", backref = 'interactor')


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
        'polymorphic_identity': 'metabolite',
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
    # remains None until determined
    is_tf = Column(String)

    __mapper_args__ = {
        'polymorphic_identity': 'protein',
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
        'polymorphic_identity': 'protein complex',
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
    strain_ortholog = Column(String)
    ortholuge_classification = Column(String)


class Interaction(Base):
    __tablename__ = 'interaction'

    id = Column(Integer, primary_key=True, autoincrement=True)
    strain = Column(String)
    # remains None unless definitively determined
    is_experimental = Column(Integer)
    # remains None unless confirmation from ortholog or derived from ortholog
    ortholog_derived = Column(String)
    comment = Column(String)
    # whether the two interactors are the same (0 or 1)
    homogenous = Column(Integer)

    xrefs = relationship("InteractionXref", backref="interaction")
    references = relationship("InteractionReference", backref="interaction")


class InteractionReference(Base):
    __tablename__ = 'interaction_reference'
    id = Column(Integer, primary_key=True, autoincrement=True)
    interaction_id = Column(Integer, ForeignKey('interaction.id'))
    pmid = Column(String)

    publication_date = Column(String)
    author_last_name = Column(String)
    publication_ref = Column(String)
    detection_method = Column(String)
    interaction_type = Column(String)
    interaction_full_name=Column(String)
    # includes type of confidence score, usually preceding confidence value with type/source:value
    confidence_score = Column(String)
    interactor_a = Column(String)
    interactor_b = Column(String)
    interactor_a_orth = Column(String)
    interactor_b_orth = Column(String)
    experimental_role_a = Column(String)
    experimental_role_b = Column(String)
    source_db = Column(String)


class InteractionXref(Base):
    __tablename__ = 'interaction_xref'

    accession = Column(String)
    interaction_id = Column(Integer, ForeignKey('interaction.id'), primary_key=True)
    data_source = Column(String, primary_key=True)






