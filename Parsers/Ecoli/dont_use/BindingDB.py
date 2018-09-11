
import csv
from Schema import Interactor, Metabolite, Interaction, InteractionReference, OrthologEcoli, InteractionXref, InteractorXref, \
    InteractionSource
from Schema import is_experimental_psimi

def parse_ecoli_bindingdb(session):
    with open('Data/Ecoli/PSICQUIC/BindingDB.txt') as csvfile:
        reader = csv.DictReader(csvfile)

        # iterate through each interaction
        for row in reader:
            uniprot_protein = None

            # check if interactor B has uniprot ID
            if 'uniprotkb' in row['ID(s) interactor B']:
                uniprot_protein = row['ID(s) interactor B'].split('uniprotkb:')[1].split('|')[0]

            if uniprot_protein is None: continue

            orthologs = []
            for ecoli_ortholog in session.query(OrthologEcoli).filter(
                    OrthologEcoli.ortholog_uniprot == uniprot_protein).all():
                if ecoli_ortholog is not None:
                    orthologs.append([ecoli_ortholog.protein, ecoli_ortholog.ortholog_id])

            if len(orthologs) == 0: continue

            ids_metabolite = row['#ID(s) interactor A'].split('|')
            chebi_metabolite, pubchem_metabolite = None, None

            # check if interactor A has ChEBI id
            for id in ids_metabolite:
                if id.split(':')[0] == 'chebi':
                    chebi_metabolite = id.split(':')[1][1:-1]

            metabolite = None

            # if interactor A has ChEBI id, query for matching metabolite
            if chebi_metabolite is not None:
                metabolite = session.query(Metabolite).filter(Metabolite.chebi == chebi_metabolite).first()

            # if unable to identify metabolite based on ChEBI id, try using pubchem id
            if metabolite is None:
                alt_ids_metabolite = row['Alt. ID(s) interactor A'].split('|')

                for id in alt_ids_metabolite:
                    if id.split(':')[0] == 'pubchem':
                        pubchem_metabolite = id.split(':')[1]

                metabolite = session.query(Metabolite).filter(Metabolite.id == pubchem_metabolite).first()

            # if unable to find interactor A in database, create new metabolite
            if metabolite is None:
                metabolite = Metabolite(id = pubchem_metabolite, pubchem = pubchem_metabolite, chebi = chebi_metabolite)
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
                    author=row['Publication 1st author(s)'].split(' ')[0]
                    date=row['Publication 1st author(s)'].split('(')[1][:-1]
                if 'pubmed:' in row['Publication Identifier(s)']:
                    pmid=row['Publication Identifier(s)'].split('pubmed:')[1][:8]

                reference = InteractionReference(interaction_id = interaction.id,
                                                 detection_method=row['Interaction detection method(s)'].split('(')[1][
                                                                  :-1],
                                                 author_ln=author,
                                                 pmid = pmid,
                                                 pub_date = date,
                                                 interaction_type = row['Interaction type(s)'].split('(')[1][:-1],
                                                 source_db=row['Source database(s)'].split('(')[1][:-1],
                                                 confidence = row['Confidence value(s)'].split('(')[0],
                                                 interactor_a=interactor_a,
                                                 interactor_b = interactor_b)

                source = session.query(InteractionSource).filter(InteractionSource.interaction_id == interaction.id,
                                                                 InteractionSource.data_source == 'BindingDB').first()

                if source is None:
                    source = InteractionSource(interaction_id=interaction.id, data_source='BindingDB')
                    session.add(source)
                session.add(reference)
        session.commit()