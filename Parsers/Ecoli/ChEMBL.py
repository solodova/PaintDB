import csv
from Schema1 import Interactor, Metabolite, Interaction, InteractionReference, OrthologEcoli, InteractionSource
from Parsers.Parser import is_experimental_psimi

def parse_ecoli_chembl(session):
    with open('Data/Ecoli/PSICQUIC/ChEMBL.txt') as csvfile:
        reader = csv.DictReader(csvfile)

        for row in reader:
            uniprot_protein = row['ID(s) interactor B'].split('uniprotkb:')[1][:6]

            orthologs = []
            for ecoli_ortholog in session.query(OrthologEcoli).filter(
                    OrthologEcoli.ortholog_uniprot == uniprot_protein).all():
                orthologs.append([ecoli_ortholog.protein, ecoli_ortholog.ortholog_id])

            if len(orthologs) == 0: continue

            chebi_metabolite = row['#ID(s) interactor A'].split('CHEBI:')[1][:6]

            metabolite = session.query(Metabolite).filter(Metabolite.chebi == chebi_metabolite).first()
            if metabolite is None:
                metabolite = Metabolite(id=chebi_metabolite)
                session.add(metabolite), session.commit()

            for interactor in orthologs:
                interaction = session.query(Interaction).filter(Interaction.interactors.contains(interactor[0]),
                                                                Interaction.interactors.contains(metabolite)).first()

                if interaction is not None:
                    if interaction.ortholog_derived is None:
                        interaction.ortholog_derived = 'cfe'
                    elif 'fe' not in interaction.ortholog_derived:
                        interaction.ortholog_derived += ', cfe'
                    session.commit()
                else:
                    interaction = Interaction(strain=interactor.strain, interactors=[metabolite, interactor[0]],
                                              type='p-m', ortholog_derived='fe')
                    # should ortholog interactions be marked as experimental?
                    if is_experimental_psimi(row['Interaction detection method(s)'].split('MI:')[1][:4]):
                        interaction.is_experimental=1
                    session.add(interaction), session.commit()

                interactor_a, interactor_b = '', ''
                if interaction.interactors[0] == metabolite:
                    interactor_a = metabolite.id
                    interactor_b = interactor[1]
                else:
                    interactor_b = metabolite.id
                    interactor_a = interactor[1]

                author, date, pmid = None, None, None
                if row['Publication 1st author(s)'] != '-':
                    author = row['Publication 1st author(s)'].split(' ')[0]
                    date = row['Publication 1st author(s)'].split('(')[1][:-1]
                if row['Publication identifier(s)'] != '-':
                    pmid = row['Publication Identifier(s)'].split('pubmed:')[1][:8]
                reference = InteractionReference(interaction_id = interaction.id,
                                                 detection_method=row['Interaction detection method(s)'].split('(')[1][
                                                                  :-1],
                                                 author_ln=author,
                                                 pmid = pmid,
                                                 pub_date = date,
                                                 interaction_type = row['Interaction type(s)'].split('(')[1][:-1],
                                                 source_db=row['Source database(s)'].split('(')[1][:-1],
                                                 confidence_score = row['Confidence value(s)'].split('(')[0],
                                                 interactor_a = interactor_a,
                                                 interactor_b = interactor_b)
                session.add(reference)
                source = session.query(InteractionSource).filter(InteractionSource.interaction_id == interaction.id,
                                                                 InteractionSource.data_source == 'ChEMBL').first()

                if source is None:
                    source = InteractionSource(interaction_id=interaction.id, data_source='ChEMBL')
                    session.add(source)

        session.commit()