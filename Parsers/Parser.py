import csv, datetime
from Parsers.PAO1_PA14 import pseudomonas_db, Parse_PSIMI_pseudomonas, regulatory_network, ortholuge_pseudomonas
from Parsers.PAO1 import Geoff_Winsor, STRING, xlinkdb, Zhang
from Parsers.PAO1_PA14_Ecoli import KEGG
from Parsers.Ecoli import EcoCyc, RegulonDB, ortholuge_ecoli, Parse_PSIMI_ecoli

# the top level function to parse all interactor info and interaction sources (ie. fill the database)
# this function is called from Main.py
def parse_all(session):
    # first obtain all info for pseudomonas proteins, then parse all interaction sources
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

    # map all existing interactions to orthologs in pseudomonas
    ortholuge_pseudomonas.parse(session)
    print(datetime.datetime.now())

    # parse Ecoli orthologs, then go through Ecoli interaction sources and map these to
    # pseudomonas orthologs
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
    # final step to update metabolite info for metabolites which may be missing ids
    KEGG.update_metabolite_info_kegg(session)
    EcoCyc.update_metabolite_info_ecocyc(session)

