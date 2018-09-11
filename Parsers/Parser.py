import csv, datetime
from Parsers.PAO1_PA14 import pseudomonas_db, Parse_PSIMI_pseudomonas, regulatory_network, ortholuge_pseudomonas
from Parsers.PAO1 import Geoff_Winsor, STRING, xlinkdb, Zhang
from Parsers.PAO1_PA14_Ecoli import KEGG
from Parsers.Ecoli import EcoCyc, RegulonDB, ortholuge_ecoli, Parse_PSIMI_ecoli

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

def parse_all(session):
    pseudomonas_db.parse(session)
    print(datetime.datetime.now())
    Geoff_Winsor.parse(session)
    print(datetime.datetime.now())
    xlinkdb.parse(session)
    print(datetime.datetime.now())
    Zhang.parse(session)
    print(datetime.datetime.now())
    regulatory_network.parse(session)
    print(datetime.datetime.now())
    Parse_PSIMI_pseudomonas.parse(session)
    print(datetime.datetime.now())
    KEGG.parse_pseudomonas(session)
    print(datetime.datetime.now())
    ortholuge_pseudomonas.parse(session)
    print(datetime.datetime.now())

    ortholuge_ecoli.parse(session)
    print(datetime.datetime.now())
    KEGG.parse_ecoli(session)
    print(datetime.datetime.now())
    EcoCyc.parse(session)
    print(datetime.datetime.now())
    Parse_PSIMI_ecoli.parse(session)
    print(datetime.datetime.now())
    RegulonDB.parse(session)
    print(datetime.datetime.now())
    KEGG.update_metabolite_info_kegg(session)
    EcoCyc.update_metabolite_info_ecocyc(session)

