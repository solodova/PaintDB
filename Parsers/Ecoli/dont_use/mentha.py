import csv
from DB_schema import OrthologEcoli, Interaction, InteractionReference, InteractionSource
from DB_build import is_experimental_psimi

def parse_ecoli_mentha(session):
    with open('Ecoli/PSICQUIC/mentha.txt') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
            if (row['#ID(s) interactor A'] == '-') | (row['ID(s) interactor B'] == '-'): continue
            interactors = []
            orthologs_A = session.query(OrthologEcoli).filter(
                OrthologEcoli.ortholog_uniprot == row['#ID(s) interactor A'].split(':')[1]).all()
            orthologs_B = session.query(OrthologEcoli).filter(
                OrthologEcoli.ortholog_uniprot == row['ID(s) interactor B'].split(':')[1]).all()

            for ortholog_A in orthologs_A:
                for ortholog_B in orthologs_B:
                    if (ortholog_A is not None) and (ortholog_B is not None):
                        if ortholog_A.strain_protein == ortholog_B.strain_protein:
                            interactors.append([[ortholog_A.protein, ortholog_A.ortholog_id],
                                                [ortholog_B.protein, ortholog_B.ortholog_id]])

            for interactor_pair in interactors:
                homogenous = (interactor_pair[0][0] == interactor_pair[1][0])
                interaction = session.query(Interaction).filter(Interaction.interactors.contains(interactor_pair[0][0]),
                                                                Interaction.interactors.contains(interactor_pair[1][0]),
                                                                Interaction.homogenous == homogenous).first()
                if interaction is not None:
                    if interaction.ortholog_derived is None:
                        interaction.ortholog_derived = 'cfe'
                    elif 'fe' not in interaction.ortholog_derived:
                        interaction.ortholog_derived += ', cfe'
                    session.commit()
                else:
                    interaction = Interaction(strain=interactor_pair[0][0].strain,
                                              interactors=[interactor_pair[0][0], interactor_pair[1][0]],
                                              type=(interactor_pair[0][0].type + '-' + interactor_pair[1][0]),
                                              ortholog_derived='fe')
                    #ask about marking ecoli ortholog interactions as experimental!!
                    if 'MI:' in row['Interaction detection method(s)']:
                        #iterate through all methods
                        if is_experimental_psimi(row['Interaction detection method(s)'].split('MI:')[1][:4]):
                            interaction.is_experimental = 1
                    session.add(interaction), session.commit()

                interactor_a, interactor_b = '', ''
                if interaction.interactors[0] == interactor_pair[0][0]:
                    interactor_a = interactor_pair[0][1]
                    interactor_b = interactor_pair[1][1]
                else:
                    interactor_b = interactor_pair[0][1]
                    interactor_a = interactor_pair[1][1]
                reference = InteractionReference(interaction_id=interaction.id,
                                                 psimi_detection=
                                                 row['Interaction detection method(s)'].split('MI:')[1][:4],
                                                 detection_method=
                                                 row['Interaction detection method(s)'].split('(')[1][:-1],
                                                 pmid=row['Publication Identifier(s)'].split('pubmed:')[1],
                                                 psimi_type=row['Interaction type(s)'].split('MI:')[1][:4],
                                                 interaction_type=row['Interaction type(s)'].split('(')[1][:-1],
                                                 psimi_db=row['Source database(s)'].split('MI:')[1][:4],
                                                 source_db = row['Source database(s)'].split('(')[1][:-1],
                                                 confidence_score=row['Confidence value(s)'])
                session.add(reference)

                source = session.query(InteractionSource).filter(
                    InteractionSource.interaction_id == interaction.id,
                    InteractionSource.data_source == 'mentha').first()

                if source is None:
                    source = InteractionSource(interaction_id=interaction.id, data_source='mentha')
                    session.add(source)

        session.commit()
        print(session.query(Interaction).count())