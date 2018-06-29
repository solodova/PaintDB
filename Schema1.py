import sqlalchemy
from sqlalchemy import create_engine, Column, Integer, String, ForeignKey, DateTime, Table
from sqlalchemy.orm import relationship, backref
from sqlalchemy.ext.declarative import declarative_base

engine = create_engine('sqlite:///:memory:', echo=True)
Base = declarative_base()

interactor_ontology = Table('interactor_ontology', Base.metadata,
                            Column('interactor_id', Integer, ForeignKey('interactor.id')),
                            Column('GO_accession', Integer, ForeignKey('ontology.GO_accession')))

interaction_participant = Table('interaction_participant', Base.metadata,
                                Column('interactor_id', Integer, ForeignKey('interactor.id')),
                                Column('interaction_id', Integer, ForeignKey('interaction.id')))

interactor_feature = Table('interactor_feature', Base.metadata,
                           Column('interactor_id', Integer, ForeignKey('interactor.id')),
                           Column('feature_id', Integer, ForeignKey('feature.id')))

interaction_reference = Table('interaction_reference', Base.metadata,
                              Column('interaction_id', Integer, ForeignKey('interaction.id')),
                              Column('reference_id', Integer, ForeignKey('reference.id')))

interaction_data_source = Table('interaction_data_source', Base.metadata,
                                Column('interaction_id', Integer, ForeignKey('interactor.id')),
                                Column('data_source_id', Integer, ForeignKey('data_source.id')))

feature_reference = Table('feature_reference', Base.metadata,
                          Column('feature_id', Integer, ForeignKey('feature.id')),
                          Column('reference_id', Integer, ForeignKey('reference.id')))


class Interactor(Base):
    # idmolecule = String/Integer?, generated/not
    # name = same as short name stored in protein/bindingsite/metabolite
    # name = Column(String)
    # # protein/metabolite/DNA binding site
    # type = Column(String, nullable=False)
    # # PAO1/PA14
    # strain = Column(String)
    # protein = relationship("Protein", uselist=False, back_populates="interactor")
    # binding_site = relationship("BindingSite", uselist=False, back_populates="binding_site")
    # metabolite = relationship("Metabolite", uselist=False, back_populates="metabolite")
    # annotation = relationship("Annotation", uselist=False, back_populates="annotation")
    # interactor_xref = relationship("InteractorXref", uselist=False, back_populates="interactor_xref")
    __tablename__ = 'interactor'
    idmolecule = Column(String, primary_key=True)
    name = Column(String)
    type = Column(String, nullable=False)
    strain = Column(String)

    __mapper_args__ = {
        'polymorphic identity': 'interactor',
        'polymorphic_on': type
    }

    xrefs = relationship("InteractorXref", backref="interactor")
    ontologies = relationship("Ontology", secondary=interactor_ontology, backref="interactors")
    interactions = relationship("Interaction", secondary=interaction_participant, backref="interactors")
    features = relationship("Feature", secondary=interactor_feature, backref="interactors")


class Protein(Interactor):
    # __tablename__ = 'protein'
    # #uniprotkb identifier?
    # id = Column(Integer, ForeignKey('interactors.idmolecule'), primary_key=True)
    # # short name for protein? maybe create an index for short name if it will be
    # # used commonly for queries
    # short_name = Column(String)
    # # full name = Column(String)
    # # strain name here?
    # interactor = relationship("Interactor", back_populates = "protein")
    __tablename__ = 'protein'
    id = Column(Integer, ForeignKey('interactors.idmolecule'), primary_key=True)
    short_name = Column(String)

    __mapper_args__ = {
        'polymorphic_identity': 'protein',
    }

    gene = relationship("Gene", back_ref=backref("protein", uselist=False))


class BindingSite(Interactor):
    # __tablename__ = 'binding_site'
    #
    # # generated id or is there binding site id?
    # id = Column(Integer, ForeignKey('interactors.idmolecule'), primary_key=True)
    # start_site = Column(Integer)
    # end_site = Column(Integer)
    # sequence = Column(String)
    # # does a binding site have a name
    # name = Column(String)
    # # strain name here?
    # interactor = relationship("Interactor", back_populates="binding_site")
    __tablename__ = 'binding_site'
    id = Column(Integer, ForeignKey('interactors.idmolecule'), primary_key=True)
    start_site = Column(Integer)
    end_site = Column(Integer)
    sequence = Column(String)
    name = Column(String)

    __mapper_args__ = {
        'polymorphic_identity': 'binding_site',
    }

    genes = relationship("Gene", backref="binding_site")


class Metabolite(Interactor):
    # __tablename__ = 'metabolite'
    # # what id should be used here
    # id = Column(Integer, primary_key=True)
    # name = Column(String)
    # # strain name here?
    # interactor = relationship("Interactor", back_populates="metabolite")
    __tablename__ = 'metabolite'
    id = Column(Integer, primary_key=True)
    name = Column(String)

    __mapper_args__ = {
        'polymorphic_identity': 'metabolite',
    }


class Gene(Base):
    __tablename__ = 'gene'
    locus_tag = Column(Integer, primary_key=True)
    name = Column(String)
    product_id = Column(Integer, ForeignKey('protein.id'))
    binding_site_id = Column(Integer, ForeignKey('binding_site.id'))

    protein = relationship("Protein", back_ref=backref("gene", uselist=False))


class Annotation(Base):
    __tablename__ = 'annotation'
    id = Column(Integer, ForeignKey('interactor.idmolecule'), primary_key=True)
    description = Column(String)

    interactor = relationship("Interactor", back_ref=backref("annotation", uselist=False))


class InteractorXref(Base):
    __tablename__ = 'interactor_xref'
    idmolecule = Column(Integer, ForeignKey("interactor.idmolecule"))
    datasource = Column(String)
    datasource_acc = Column(String, primary_key=True)


class Ontology(Base):
    __tablename__ = 'ontology'
    GO_accession = Column(Integer, primary_key=True)
    description = Column(String)


class Feature(Base):
    __tablename__ = 'feature'
    id = Column(Integer, primary_key=True)
    name = Column(String)


class Reference(Base):
    __tablename__ = 'reference'
    id = Column(Integer, primary_key=True)
    publication_date = Column(DateTime)
    author_first_name = Column(String)


class DataSource(Base):
    __tablename__ = 'data_source'
    id = Column(String, primary_key=True)
    name = Column(String)


class Interaction(Base):
    __tablename__ = 'interaction'
    id = Column(Integer, primary_key=True)
    type = Column(String)
    specific_type = Column(String)
    strain = Column(String)
    evidence = Column(String)

    xrefs = relationship("InteractionXref", back_ref="interaction")
    data_sources = relationship("DataSource", secondary=interaction_data_source, backref="interactions")
    references = relationship("Reference", secondary=interaction_reference, backref="interactions")


class InteractionXref(Base):
    __tablename__ = 'interaction_xref'
    interaction_id = Column(Integer, ForeignKey("interaction.id"))
    data_source_id = Column(String)
    data_source_accession = Column(String, primary_key=True)
