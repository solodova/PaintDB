import sqlalchemy
from sqlalchemy import create_engine, Column, Integer, String, ForeignKey, DateTime, Table, Float
from sqlalchemy.orm import relationship
from sqlalchemy.orm import backref
from sqlalchemy.ext.declarative import declarative_base


Base = declarative_base()

interactor_reference = Table('interactor_reference', Base.metadata,
                             Column('interactor_id', String, ForeignKey('interactor.id')),
                             Column('reference_pmid', String, ForeignKey('reference.pmid')))

interaction_participant = Table('interaction_participant', Base.metadata,
                                Column('interactor_id', String, ForeignKey('interactor.id')),
                                Column('interaction_id', Integer, ForeignKey('interaction.id')))
# feature_reference = Table('feature_reference', Base.metadata,
#                           Column('feature_id', Integer, ForeignKey('feature.id')),
#                           Column('reference_id', Integer, ForeignKey('reference.PMID')))
#
# database_reference = Table('database_reference', Base.metadata,
#                               Column('database', String, ForeignKey('database.id')),
#                               Column('interaction_reference', Integer, ForeignKey('interaction_reference.id')))


class Interactor(Base):
    __tablename__ = 'interactor'
    id = Column(String, primary_key=True)

    name = Column(String)
    type = Column(String)
    strain = Column(String)
    description = Column(String)
    accession = Column(String)
    uniprotkb = Column(String)
    start_site = Column(String)
    end_site = Column(String)
    length = Column(String)
    #just for transcription factors
    binding_sequence = Column(String)

    xrefs = relationship("InteractorXref", backref="interactor")
    gene_ontologies = relationship("GeneOntology", backref="interactor")
    interactions = relationship("Interaction", secondary=interaction_participant, backref="interactors")
    #features = relationship("Feature", backref="interactor")
    references = relationship("Reference", secondary=interactor_reference, backref="interactors")
    #binding_site_genes = relationship("Interactor")
    localizations = relationship("Localization", backref="interactor")
    orthologs = relationship("Ortholog", backref = "interactor")

class InteractorXref(Base):
    __tablename__ = 'interactor_xref'
    accession = Column(String, primary_key=True)
    interactor_id = Column(String, ForeignKey("interactor.id"), primary_key=True)

    data_source = Column(String)


class GeneOntology(Base):
    __tablename__ = 'gene_ontology'
    id = Column(Integer, primary_key=True, autoincrement=True)
    interactor_id = Column(String, ForeignKey('interactor.id'))

    accession = Column(String)
    ontology = Column(String)
    go_term = Column(String)
    evidence_code = Column(String)
    evidence_description = Column(String)
    eco_code = Column(String)
    eco_term = Column(String)
    pmid = Column(String)


class Reference(Base):
    __tablename__ = 'reference'
    pmid = Column(String, primary_key=True)

    publication_date = Column(String)
    author_last_name = Column(String)
    publication_ref = Column(String)


class Localization(Base):
    __tablename__ = 'localization'
    id = Column(Integer, primary_key=True, autoincrement=True)
    interactor_id = Column(String, ForeignKey('interactor.id'))

    localization = Column(String)
    confidence = Column(String)
    pmid = Column(String)


class Interaction(Base):
    __tablename__ = 'interaction'
    id = Column(Integer, primary_key=True, autoincrement=True)

    type = Column(String)
    strain = Column(String)
    is_experimental = Column(Integer)
    is_Ecoli_ortholog = Column(Integer)
    #num_interactors = Column(Integer)

    xrefs = relationship("InteractionXref", backref="interaction")
    references = relationship("InteractionReference", backref="interactions")
    #features = relationship("Feature", backref="interaction")


# class Feature(Base):
#     __tablename__ = 'feature'
#     id = Column(String, primary_key=True)
#     interaction_id = Column(Integer, ForeignKey('interaction.id'))
#     interactor_id = Column(String, ForeignKey('interactor.id'))
#
#     interactor_name = Column(String)
#     name = Column(String)
#     reference_PMID = Column(String)


# class Database(Base):
#     __tablename__ = 'database'
#     id = Column(String, primary_key=True)
#
#     name = Column(String)
#
#     references = relationship("InteractionReference",secondary=database_reference, backref="databases")
#

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
    accession = Column(String, primary_key=True)
    interaction_id = Column(Integer, ForeignKey('interaction.id'))
    data_source = Column(String)

class OrthologPseudomonas(Base):
    __tablename__ = 'ortholog_pseudomonas'

    interactor_id  = Column(String, ForeignKey('interactor.id'), primary_key=True)
    ortholog_id = Column(String, primary_key=True)
    strain_interactor = Column(String)
    strain_ortholog = Column(String)

class OrthologEcoli(Base):
    __tablename__ = 'ortholog_ecoli'

    interactor_id = Column(String, ForeignKey('interactor.id'), primary_key=True)
    ortholog_id = Column(String, primary_key=True)
    ortholog_uniprot = Column(String)
    ortholog_refseq = Column(String)
    ortholog_name = Column(String)
    strain_interactor = Column(String)
    strain_ortholog = Column(String)
    ortholuge_classification = Column(String)





