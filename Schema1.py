import sqlalchemy
from sqlalchemy import create_engine, Column, Integer, String, ForeignKey, DateTime, Table, Float
from sqlalchemy.orm import relationship
from sqlalchemy.orm import backref
from sqlalchemy.ext.declarative import declarative_base

experimental_psimi = ['0045', '0401', '0013', '0254', '0428', '1088', '0255', '0090', '0400', '1232', '0091', '0027',
                      '0030', '0982', '0415', '0417', '1036', '0953', '2197', '0968', '0016', '0943', '1311', '0894',
                      '0043', '1219', '1086', '0928', '0051', '0964', '0859', '0065', '0067', '2169', '1247', '0071',
                      '0893', '0891', '0077', '0938', '0099', '0888', '1235', '0966', '0114', '2224', '0439', '0441',
                      '0872', '0663', '0040', '0416', '0426', '2213', '0827', '1089', '1192', '0257', '2285', '0256',
                      '0014', '0010', '0011', '0809', '0111', '0231', '0895', '0097', '1203', '0370', '0232', '0004',
                      '0008', '0405', '0034', '1087', '1031', '0813', '0440', '0892', '0657', '0226', '0227', '0028',
                      '0029', '1022', '1211', '0430', '0031', '0807', '0983', '0889', '1005', '0989', '1142', '1147',
                      '1137', '0990', '1309', '0406', '2281', '0984', '2216', '1138', '0996', '1006', '0870', '1011',
                      '1009', '0998', '1026', '0999', '1007', '1000', '1249', '1001', '0515', '1010', '1034', '0879',
                      '0979', '0434', '1145', '0972', '0841', '0696', '1325', '1008', '0997', '1229', '0602', '0605',
                      '1191', '0949', '2198', '0969', '1342', '2196', '1038', '0107', '0069', '0944', '0041', '0042',
                      '0905', '0229', '0012', '0017', '1030', '0052', '0053', '1016', '0054', '0055', '0510', '0976',
                      '0965', '0038', '0104', '0826', '2170', '2171', '1104', '1103', '0425', '0825', '0824', '0410',
                      '1024', '0020', '1037', '1204', '0655', '0369', '1320', '0432', '0726', '0588', '0018', '2288',
                      '0019', '0096', '0676', '0225', '0081', '0089', '0921', '0073', '0084', '0098', '0115', '2189',
                      '0947', '0411', '0047', '0049', '2167', '0729', '2191', '2194', '1312', '0404', '0808', '0413',
                      '0887', '1354', '0991', '0435', '0508', '1252', '1236', '1003', '1002', '1004', '0516', '1035',
                      '0920', '0880', '0419', '0514', '1019', '0424', '0697', '0698', '0699', '0700', '0603', '1190',
                      '0604', '0901', '1189', '1183', '0814', '1352', '1313', '2199', '1246', '1238', '0009', '0420',
                      '0509', '0511', '1321', '0112', '0437', '0438', '0728', '0727', '0916', '0397', '0398', '1356',
                      '0006', '0007', '0402', '0858', '0946', '1017', '0963', '0678', '0092', '0095', '0048', '0066',
                      '0108', '2283', '2192', '2193', '2188', '2195', '0276', '0412', '0992', '0993', '0994', '0995',
                      '0513', '0512', '0423', '0606', '2168', '1184', '1314', '1113', '1111', '1218', '1028', '1029',
                      '0695', '0899', '0900', '1187', '2277', '2215', '1112', '0399']


def is_experimental_psimi(psi_code):
    return psi_code in experimental_psimi

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

    localizations = relationship("Localization", secondary=protein_localizations, backref="protein")
    pseudomonas_orthologs = relationship("OrthologPseudomonas", backref ="protein")
    ecoli_orthologs = relationship("OrthologEcoli", backref = "protein")
    xrefs = relationship("ProteinXref", backref='protein')


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

    id = Column(Integer, primary_key=True, autoincrement=True)
    data_source=Column(String)
    is_experimental = Column(Integer)






