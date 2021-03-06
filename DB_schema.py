from sqlalchemy import Column, Integer, String, ForeignKey, Table
from sqlalchemy.orm import relationship
from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()

protein_localizations = Table('protein_localizations', Base.metadata,
                              Column('localization_id', String, ForeignKey('localization.id')),
                              Column('protein_id', String, ForeignKey('protein.id')))

interaction_participants = Table('interaction_participants', Base.metadata,
                                Column('interactor_id', String, ForeignKey('interactor.id')),
                                Column('interaction_id', Integer, ForeignKey('interaction.id')))

interaction_references = Table('interaction_references', Base.metadata,
                                Column('interaction_id', String, ForeignKey('interaction.id')),
                                Column('reference_id', Integer, ForeignKey('interaction_reference.id')))

interaction_sources = Table('interaction_sources', Base.metadata,
                                Column('interaction_id', String, ForeignKey('interaction.id')),
                                Column('data_source', Integer, ForeignKey('interaction_source.id')))

reference_sources = Table('reference_sources', Base.metadata,
                          Column('reference_id', String, ForeignKey('interaction_reference.id')),
                          Column('data_source', String, ForeignKey('interaction_source.id')))


class Interactor(Base):
    __tablename__ = 'interactor'

    # protein: locus id (eg. PA0001)
    # protein complex: uniprotkb id (eg. Q9HT76)
    # metabolite: any one of the ids by which it was defined when parsing (can be kegg, pubchem, ecocyc, cas, chebi)
    # Note: search for metabolites by  specific attributes rather than id
    id = Column(String, primary_key=True)
    # protein: eg. pyrE (None until assigned)
    # protein complex: None
    # metabolite: name which was found while parsing (None until assigned)
    name = Column(String)
    type = Column(String)

    __mapper_args__ = {
        'polymorphic_identity': 'interactor',
        'polymorphic_on': type
    }

    interactions = relationship("Interaction", secondary=interaction_participants, backref='interactors')


class Metabolite(Interactor):
    __tablename__ = 'metabolite'

    #can be any one of the ids listed below
    # Note: query metabolite by specific attribute rather than id because of this inconsistency
    # if the database is rerun, can change this to an autoincremented value
    id = Column(String, ForeignKey('interactor.id'), primary_key=True)
    kegg = Column(String)
    pubchem = Column(String)
    cas = Column(String)
    chebi = Column(String)
    ecocyc = Column(String)

    __mapper_args__ = {
        'polymorphic_identity': 'm',
    }


class Protein(Interactor):
    __tablename__ = 'protein'

    # monomeric protein: locus id (eg. PA0001)
    # protein complex: uniprot id
    id = Column(String, ForeignKey('interactor.id'), primary_key=True)
    strain = Column(String)
    # monomeric protein: product name (eg. orotate phosphoribosyltransferase, hypothetical protein)
    # protein complex: None
    product_name = Column(String)
    # eg. NP_254184.1 (None until specified)
    ncbi_acc = Column(String)
    # eg. Q9HT76 (None until specified)
    uniprotkb = Column(String)

    __mapper_args__ = {
        'polymorphic_identity': 'p',
    }

    localizations = relationship("Localization", secondary=protein_localizations, backref="proteins")
    pseudomonas_orthologs = relationship("OrthologPseudomonas", backref ="protein")
    ecoli_orthologs = relationship("OrthologEcoli", backref = "protein")
    ontologies = relationship("GeneOntology", backref = 'protein')


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
    # eg. Cytoplasmic
    localization = Column(String)
    # eg. Class 3
    confidence = Column(String)


class OrthologPseudomonas(Base):
    __tablename__ = 'ortholog_pseudomonas'

    # id of the protein interactor
    # Note: a protein may have >1 ortholog in the other pseudomonas strain
    protein_id  = Column(String, ForeignKey('protein.id'), primary_key=True)
    # id of the ortholog in the opposite strain
    ortholog_id = Column(String, primary_key=True)
    # strain of the protein interactor
    strain_protein = Column(String, primary_key=True)
    # strain of the ortholog interactor
    strain_ortholog = Column(String)


class OrthologEcoli(Base):
    __tablename__ = 'ortholog_ecoli'

    # id of the protein interactor
    # Note: a protein may have >1 ortholog in the Ecoli
    protein_id = Column(String, ForeignKey('protein.id'), primary_key=True)
    # strain of protein interactor
    strain_protein = Column(String)
    # ecoli ortholog info
    ortholog_id = Column(String, primary_key=True)
    ortholog_uniprot = Column(String)
    ortholog_refseq = Column(String)
    ortholog_name = Column(String)
    # classification given by ortholuge analysis (will be RBBH or SSD)
    ortholuge_classification = Column(String)


class Interaction(Base):
    __tablename__ = 'interaction'

    id = Column(Integer, primary_key=True, autoincrement=True)
    strain = Column(String)
    # p-p, m-p, or p-m
    type = Column(String)
    # whether the two interactors are the same (0 or 1)
    homogenous = Column(Integer)
    # strain the interaction was derived from (if not found in source for its own strain)
    ortholog_derived = Column(String)

    xrefs = relationship("InteractionXref", backref="interaction")
    references = relationship("InteractionReference", secondary=interaction_references, backref="interactions")
    sources = relationship("InteractionSource", secondary=interaction_sources, backref="interactions")


class InteractionReference(Base):
    __tablename__ = 'interaction_reference'
    id = Column(Integer, primary_key=True, autoincrement=True)
    #interaction_id = Column(Integer, ForeignKey('interaction.id'))

    # all values here are None until specified

    # may be PSIMI term, may not be
    detection_method = Column(String)
    author_ln = Column(String)
    pub_date = Column(String)
    pmid = Column(String)
    # may be PSIMI or not
    interaction_type = Column(String)
    # may be PSIMI or not
    source_db = Column(String)
    # when there is additional comment about the interaction (eg. more detailed description)
    comment=Column(String)
    # includes type of confidence score, usually preceding confidence value with type/source:value
    confidence = Column(String)
    # interactor a and b ids if the reference is from an ortholog source (specifies the interactors
    # in the original interaction)
    interactor_a = Column(String)
    interactor_b = Column(String)

    sources = relationship('InteractionSource', secondary=reference_sources, backref='references')


class InteractionXref(Base):
    __tablename__ = 'interaction_xref'

    accession = Column(String, primary_key=True)
    interaction_id = Column(Integer, ForeignKey('interaction.id'), primary_key=True)
    data_source = Column(String)


#source is considered to be one of the sources from which information was parsed (needed for query filter)
class InteractionSource(Base):
    __tablename__ ='interaction_source'

    id = Column(Integer, primary_key=True, autoincrement=True)
    data_source=Column(String)
    is_experimental = Column(Integer)






