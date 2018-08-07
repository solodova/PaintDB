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

    # locus id (protein)
    # uniprotkb id (protein complex)
    # kegg id (metabolite
    id = Column(String, primary_key=True)
    name = Column(String)
    type = Column(String)

    __mapper_args__ = {
        'polymorphic_identity': 'interactor',
        'polymorphic_on': type
    }

    interactions = relationship("Interaction", secondary=interaction_participant, backref="interactors")


class Metabolite(Interactor):
    __tablename__ = 'metabolite'

    __mapper_args__ = {
        'polymorphic_identity': 'metabolite',
    }

    #KEGG, then EcoCyc code
    id = Column(String, ForeignKey('interactor.id'), primary_key=True)
    KEGG = Column(String)
    PubChem = Column(String)
    CAS = Column(String)
    ChEBI = Column(String)
    EcoCyc = Column(String)


class Protein(Interactor):
    __tablename__ = 'protein'

    # locus id (monomeric protein) or uniprotkb id (protein complex)
    id = Column(String, ForeignKey('interactor.id'), primary_key=True)
    # description of protein product for monomeric, 'protein complex' otherwise
    description = Column(String)
    is_TF = Column(String)
    accession = Column(String)
    uniprotkb = Column(String)
    strain = Column(String)

    __mapper_args__ = {
        'polymorphic_identity': 'protein',
    }

    localizations = relationship("Localization", backref="protein")
    references = relationship("Reference", secondary=protein_reference, backref="proteins")
    xrefs = relationship("ProteinXref", backref="protein")
    pseudomonas_orthologs = relationship("OrthologPseudomonas", backref ="protein")
    ecoli_ortholgs = relationship("OrthologEcoli", backref = "protein")


class ProteinXref(Base):
    __tablename__ = 'protein_xref'

    accession = Column(String, primary_key=True)
    protein_id = Column(String, ForeignKey("protein.id"), primary_key=True)
    data_source = Column(String)


class Reference(Base):
    __tablename__ = 'reference'

    pmid = Column(String, primary_key=True)
    publication_date = Column(String)
    author_last_name = Column(String)
    publication_ref = Column(String)


class GeneOntology(Base):
    __tablename__ = 'gene_ontology'

    # automatically generated
    id = Column(Integer, primary_key=True, autoincrement=True)
    gene_id = Column(String, ForeignKey('protein.id'))
    # GO accession
    accession = Column(String)
    go_term = Column(String)
    evidence_code = Column(String)
    eco_code = Column(String)
    eco_term = Column(String)
    pmid = Column(String)


class Localization(Base):
    __tablename__ = 'localization'

    id = Column(Integer, primary_key=True, autoincrement=True)
    protein_id = Column(String, ForeignKey('protein.id'))
    localization = Column(String)
    confidence = Column(String)
    pmid = Column(String)


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
    type = Column(String)
    strain = Column(String)
    is_experimental = Column(Integer)
    ortholog_derived = Column(String)
    comment = Column(String)
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
    interactor_a_id = Column(String)
    interactor_b_id = Column(String)
    experimental_role_a = Column(String)
    experimental_role_b = Column(String)
    source_db = Column(String)


class InteractionXref(Base):
    __tablename__ = 'interaction_xref'

    accession = Column(String)
    interaction_id = Column(Integer, ForeignKey('interaction.id'), primary_key=True)
    data_source = Column(String, primary_key=True)






