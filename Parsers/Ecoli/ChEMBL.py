import csv
from Schema1 import Interactor, Metabolite, Interaction, InteractionReference, OrthologEcoli, InteractionXref

def parse_Ecoli_ChEMBL(session):
    with open('Data/Ecoli/PSICQUIC/ChEMBL.txt') as csvfile:
        reader = csv.DictReader(csvfile)

        for row in reader:
            IDs_metabolite = row['#ID(s) interactor A'].split('|')
            IDs_protein = row['ID(s) interactor B'].split('|')
            ChEBI_metabolite = None
            UniProt_protein = None
            for id in IDs_metabolite:
                if id[:5] == 'chebi':
                    ChEBI_metabolite = id.split(':')[2][:-1]

            for id in IDs_protein:
                if id[:8] == 'uniprotkb':
                    UniProt_protein = id.split(':')[1]


            metabolite = session.query(Metabolite).filter(Metabolite.ChEBI == ChEBI_metabolite).first()
            if metabolite is None: continue

            interactor_B = []
            for ortholog in session.query(OrthologEcoli).filter(OrthologEcoli.ortholog_uniprot == UniProt_protein).all():
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