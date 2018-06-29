import sqlalchemy
from sqlalchemy import create_engine, Column, Integer, String, ForeignKey, DateTime, Table, Float
from sqlalchemy.orm import relationship, backref
from sqlalchemy.ext.declarative import declarative_base

engine = create_engine('sqlite:///:memory:', echo=True)
Base = declarative_base()

interactor_reference = Table('interactor_reference', Base.metadata,
                             Column('interactor_id', String, ForeignKey('interactor.id')),
                             Column('reference_id'), Integer, ForeignKey('reference.PMID'))

interaction_participant = Table('interaction_participant', Base.metadata,
                                Column('interactor_id', String, ForeignKey('interactor.id')),
                                Column('interaction_id', Integer, ForeignKey('interaction.id')))

feature_reference = Table('feature_reference', Base.metadata,
                          Column('feature_id', Integer, ForeignKey('feature.id')),
                          Column('reference_id', Integer, ForeignKey('reference.PMID')))

data_source_reference = Table('data_source_reference', Base.metadata,
                              Column('data_source', String, ForeignKey('data_source.id')),
                              Column('interaction_reference', Integer, ForeignKey('interaction_reference.id')))


class Interactor(Base):
    __tablename__ = 'interactor'
    id = Column(String, primary_key=True)

    name = Column(String)
    type = Column(String, nullable=False)
    strain = Column(String)
    description = Column(String)

    __mapper_args__ = {
        'polymorphic identity': 'interactor',
        'polymorphic_on': type
    }

    xrefs = relationship("InteractorXref", backref="interactor")
    gene_ontologies = relationship("GeneOntology", backref="interactor")
    interactions = relationship("Interaction", secondary=interaction_participant, backref="interactors")
    features = relationship("Feature", backref="interactor")
    references = relationship("Reference", secondary=interactor_reference, backref="interactors")


class BindingSite(Interactor):
    __tablename__ = 'binding_site'
    id = Column(Integer, ForeignKey('interactors.id'), primary_key=True)

    start_site = Column(Integer)
    end_site = Column(Integer)
    sequence = Column(String)

    __mapper_args__ = {
        'polymorphic_identity': 'binding_site',
    }

    genes = relationship("Gene", backref="binding_site")


class Gene(Base):
    __tablename__ = 'gene'
    locus_tag = Column(Integer, primary_key=True)
    protein_id = Column(Integer, ForeignKey('interactor.id'))
    binding_site_id = Column(Integer, ForeignKey('binding_site.id'))

    name = Column(String)
    start_site = Column(Integer)
    end_site = Column(Integer)

    protein = relationship("Interactor", back_ref=backref("gene", uselist=False))


class InteractorXref(Base):
    __tablename__ = 'interactor_xref'
    accession = Column(String, primary_key=True)
    interactor_id = Column(Integer, ForeignKey("interactor.id"))

    data_source = Column(String)


class GeneOntology(Base):
    __tablename__ = 'gene_ontology'
    id = Column(String, primary_key=True)
    interactor_id = Column(String, ForeignKey('interactor.id'))

    accession = Column(String)
    ontology = Column(String)
    term = Column(String)
    evidence_code = Column(String)
    evidence_description = Column(String)
    ECO_code = Column(String)
    ECO_description = Column(String)
    PMID = Column(Integer)


class Reference(Base):
    __tablename__ = 'reference'
    PMID = Column(Integer, primary_key=True)

    publication_date = Column(DateTime)
    author_last_name = Column(String)


class Localization(Base):
    __tablename__ = 'localization'
    interactor_id = Column(String, ForeignKey('interactor.id'), primary_key=True)

    localization = Column(String)
    confidence = Column(String)
    PMID = Column(Integer)

    interactor = relationship("Interactor", back_ref=backref("localization", uselist=False))


class Interaction(Base):
    __tablename__ = 'interaction'
    id = Column(Integer, primary_key=True)

    type = Column(String)
    strain = Column(String)

    xrefs = relationship("InteractionXref", back_ref="interaction")
    data_sources = relationship("DataSource", backref="interactions")
    references = relationship("InteractionReference", backref="interactions")
    features = relationship("Feature", backref="interaction")

class Feature(Base):
    __tablename__ = 'feature'
    id = Column(Integer, primary_key=True)
    interaction_id = Column(Integer, ForeignKey('interaction.id'))
    interactor_id = Column(String, ForeignKey('interactor.id'))

    interactor_name = Column(String)
    name = Column(String)
    reference_PMID = Column(Integer)


class DataSource(Base):
    __tablename__ = 'data_source'
    id = Column(String, primary_key=True)
    interaction_id = Column(Integer, ForeignKey('interaction.id'))

    db_name = Column(String)
    accession = Column(String)

    references = relationship("InteractionReference",secondary=data_source_reference, backref="data_sources")


class InteractionReference(Base):
    __tablename__ = 'interaction_reference'
    id = Column(Integer, primary_key=True)
    interaction_id = Column(Integer, ForeignKey('interaction.id'))

    PMID = Column(Integer)
    publication_date = Column(DateTime)
    author_last_name = Column(String)
    detection_method = Column(String)
    interaction_type = Column(String)
    confidence_score = Column(Float)
    interactor_A = Column(String)
    interactor_B = Column(String)
    interactor_A_id = Column(String)
    interactor_B_id = Column(String)
    experimental_role_A = Column(String)
    experimental_role_B = Column(String)


class InteractionXref(Base):
    __tablename__ = 'interaction_xref'
    accession = Column(String, primary_key=True)
    interaction_id = Column(Integer, ForeignKey("interaction.id"))

    data_source = Column(String)





