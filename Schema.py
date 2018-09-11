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

    interactions = relationship("Interaction", secondary=interaction_participants, backref='interactors')


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

    localizations = relationship("Localization", secondary=protein_localizations, backref="proteins")
    pseudomonas_orthologs = relationship("OrthologPseudomonas", backref ="protein")
    ecoli_orthologs = relationship("OrthologEcoli", backref = "protein")
    xrefs = relationship("ProteinXref", backref='protein')
    ontologies = relationship("GeneOntology", backref = 'protein')


class ProteinXref(Base):
    # note that some PAO1 proteins will have multiple xrefs from the same source
    __tablename__ = 'interactor_xref'

    protein_id = Column(String, ForeignKey("interactor.id"), primary_key=True)
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

    protein_id  = Column(String, ForeignKey('protein.id'), primary_key=True)
    ortholog_id = Column(String, primary_key=True)
    strain_protein = Column(String, primary_key=True)
    strain_ortholog = Column(String)


class OrthologEcoli(Base):
    __tablename__ = 'ortholog_ecoli'

    protein_id = Column(String, ForeignKey('protein.id'), primary_key=True)
    strain_protein = Column(String)
    ortholog_id = Column(String, primary_key=True)
    ortholog_uniprot = Column(String)
    ortholog_refseq = Column(String)
    ortholog_name = Column(String)
    ortholuge_classification = Column(String)


class Interaction(Base):
    __tablename__ = 'interaction'

    id = Column(Integer, primary_key=True, autoincrement=True)
    strain = Column(String)
    type = Column(String)
    # whether the two interactors are the same (0 or 1)
    homogenous = Column(Integer)
    ortholog_derived = Column(String)

    xrefs = relationship("InteractionXref", backref="interaction")
    references = relationship("InteractionReference", secondary=interaction_references, backref="interactions")
    sources = relationship("InteractionSource", secondary=interaction_sources, backref="interactions")


class InteractionReference(Base):
    __tablename__ = 'interaction_reference'
    id = Column(Integer, primary_key=True, autoincrement=True)
    #interaction_id = Column(Integer, ForeignKey('interaction.id'))

    # PSI MI term (except Geoff)
    detection_method = Column(String)
    author_ln = Column(String)
    pub_date = Column(String)
    pmid = Column(String)
    # remains None if '-', otherwise psimi OR just 'predicted'
    interaction_type = Column(String)
    source_db = Column(String)
    # when there is additional comment about the interaction (eg. more detailed description)
    comment=Column(String)
    # includes type of confidence score, usually preceding confidence value with type/source:value
    confidence = Column(String)
    interactor_a = Column(String)
    interactor_b = Column(String)

    sources = relationship('InteractionSource', secondary=reference_sources, backref='references')


class InteractionXref(Base):
    __tablename__ = 'interaction_xref'

    accession = Column(String, primary_key=True)
    interaction_id = Column(Integer, ForeignKey('interaction.id'), primary_key=True)
    data_source = Column(String)


class InteractionSource(Base):
    __tablename__ ='interaction_source'

    id = Column(Integer, primary_key=True, autoincrement=True)
    data_source=Column(String)
    is_experimental = Column(Integer)






