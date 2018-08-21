
import csv
from Schema1 import Interactor, Metabolite, Interaction, InteractionReference, OrthologEcoli, InteractionXref, InteractorXref

def parse_Ecoli_BindingDB(session):
    with open('Data/Ecoli/PSICQUIC/BindingDB.txt') as csvfile:
        reader = csv.DictReader(csvfile)

        # iterate through each interaction
        for row in reader:
            ids_A = row['#ID(s) interactor A'].split('|')
            ChEBI_A = None
            pubchem_A = None
            BindingDB_A = None

            # check if interactor A has ChEBI id
            for id in ids_A:
                if id.split(':')[0] == 'chebi':
                    ChEBI_A = id.split(':')[1][1:-1]
                elif id.split(':')[0] == 'BindingDB_monomerID':
                    BindingDB_A = id.split(':')[1]

            interactor_A = None

            # if interactor A has ChEBI id, query for matching metabolite
            if ChEBI_A is not None:
                interactor_A = session.query(Metabolite).filter(Metabolite.ChEBI == ChEBI_A).first()

            # if unable to identify metabolite based on ChEBI id, try using pubchem id
            if interactor_A is None:
                alt_IDs_A = row['Alt. ID(s) interactor A'].split('|')

                for id in alt_IDs_A:
                    if id.split(':')[0] == 'pubchem':
                        pubchem_A = id.split(':')[1]

                interactor_A = session.query(Metabolite).filter(Metabolite.id == pubchem_A).first()

            # if unable to find interactor A in database, create new metabolite
            if (interactor_A == None):
                interactor_A = Metabolite(id = pubchem_A, PubChem = pubchem_A, ChEBI = ChEBI_A)
                #create bindingdb xref
                interactor_A.xrefs.append(InteractorXref(interactor_id = pubchem_A, accession =BindingDB_A,
                                                         source = 'BindingDB_monomerID'))
                session.add(interactor_A)
                session.commit()

            IDs_B = row['ID(s) interactor B'].split('|')
            interactor_B = None
            UniProt_B = None
            # check if interactor B has uniprot ID
            for id in IDs_B:
                if id[:5] == 'uniprotkb':
                    UniProt_B = id.split(':')[1]

            # if interactor B has uniprot id, try to find it in database
            if UniProt_B is not None:
                interactor_B = session.query(OrthologEcoli).filter(OrthologEcoli.ortholog_uniprot == UniProt_B).first()

            # if could not find interactor B in database, try alt. ids (also uniprot)
            if interactor_B is None:
                for id in row['Alt. ID(s) interactor B']:
                    if id[:8] == 'uniprotkb':
                        UniProt_B = id.split(':')[1]

                # if interactor has uniprot ID in alt ids, try to find it in database
                if UniProt_B is not None:
                    interactor_B = session.query(OrthologEcoli).filter(OrthologEcoli.ortholog_uniprot == UniProt_B).first()

            # if interactor B was not found, go to next interaction, since no orthologs for this protein exist
            if interactor_B is None: continue

            orthologs_B = []
            for ecoli_protein in session.query(OrthologEcoli).filter(OrthologEcoli.ortholog_uniprot == UniProt_protein).all():
                interactor_B.append([ortholog.protein, ortholog.ortholog_id])

            for interactor in interactor_B:
                interaction = session.query(Interaction).filter(Interaction.interactors.contains(interactor[0]),
                                                                Interaction.interactors.contains(metabolite)).first()

                if (interaction != None):
                    if (interaction.ortholog_derived == None):
                        interaction.ortholog_derived = 'confirmed from E.coli'
                    elif ('from E. coli' not in interaction.ortholog_derived):
                        interaction.ortholog_derived += ', confirmed from E. coli'
                    session.commit()
                else:
                    interaction = Interaction(strain=interactor.strain, interactors=[metabolite, interactor[0]],
                                              type='protein-metabolite',
                                              is_experimental=0, ortholog_derived='from E. coli')
                    session.add(interaction), session.commit()

                reference = InteractionReference(interaction_id = interaction.id,
                                                 pmid = row['Publication Identifier(s)'].split(':')[1],
                                                 publication_date = row['Publication 1st author(s)'].split('(')[1][:-1],
                                                 author_last_name = row['Publication 1st author(s)'].split(' ')[0],
                                                 publication_ref = row['Publication 1st author(s)'],
                                                 detection_method = row['Interaction detection method(s)'].split('(')[1][:-1],
                                                 interaction_type = row['Interaction type(s)'].split('(')[1][:-1],
                                                 confidence_score = row['Confidence value(s)'].split('(')[0],
                                                 interactor_A = metabolite.id, interactor_B = interactor.id,
                                                 interactor_A_id = metabolite.id, interactor_B_id = interactor[1],
                                                 source_db = row['Source database(s)'].split('(')[1][:-1])
                xref = InteractionXref(interaction_id = interaction.id, data_source = 'ChEMBL',
                                       accession = row['Interaction identifier(s)'].split(':')[1])
                session.add(reference), session.add(xref), session.commit()