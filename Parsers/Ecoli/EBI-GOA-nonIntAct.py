import csv
from Schema1 import Interactor, Interaction, OrthologEcoli, InteractionReference, InteractionSource
from Parsers.Parser import is_experimental_psimi

def parse_ecoli_ebi_goa_nonintact(session):
    with open('Ecoli/PSICQUIC/EBI-GOA-nonIntAct.txt') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
            interactors = []
            uniprot_A, uniprot_B = None, None
            if 'uniprotkb:' in row['#ID(s) interactor A']:
                uniprot_A = row['#ID(s) interactor A'].split('uniprotkb:')[1]
            if 'uniprotkb:' in row['ID(s) interactor B']:
                uniprot_B = row['ID(s) interactor B'].split('uniprotkb:')[1]

            if (uniprot_A is None) | (uniprot_B is None): continue

            orthologs_A = session.query(OrthologEcoli).filter(OrthologEcoli.ortholog_uniprot == uniprot_A).all()
            orthologs_B = session.query(OrthologEcoli).filter(OrthologEcoli.ortholog_uniprot == uniprot_B).all()
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
                    interaction = Interaction(strain=interactor_pair[0][0].strain, interactors=interactor_pair,
                                              type='p-p', ortholog_derived='fe')
                    if is_experimental_psimi(row['Interaction detection method(s)'].split('MI:')[1][:4]):
                        interaction.is_experimental=1
                    session.add(interaction), session.commit()

                interactor_a, interactor_b = '', ''
                if interaction.interactors[0] == interactor_pair[0][0]:
                    interactor_a = interactor_pair[0][1]
                    interactor_b = interactor_pair[1][1]
                else:
                    interactor_b = interactor_pair[0][1]
                    interactor_a = interactor_pair[1][1]

                reference = InteractionReference(interaction_id=interaction.id,
                                                 detection_method=row['Interaction detection method(s)'].split('(')[
                                                                      1][:-1],
                                                 author_ln=row['Publication 1st author(s)'].split(' ')[0],
                                                 pub_date=row['Publication 1st author(s)'].split('(')[1],
                                                 pmid=row['Publication Identifier(s)'].split('pubmed:')[1],
                                                 interaction_type=row['Interaction type(s)'].split('(')[1][:-1],
                                                 source_db=row['Source database(s)'].split('(')[1][:-1],
                                                 interactor_a_id=row['#ID(s) interactor A'].split(':')[1],
                                                 interactor_b_id=row['ID(s) interactor B'].split(':')[1])
                session.add(reference)

                source = session.query(InteractionSource).filter(
                    InteractionSource.interaction_id == interaction.id,
                    InteractionSource.data_source == 'EBI-GOA non-IntAct').first()

                if source is None:
                    source = InteractionSource(interaction_id=interaction.id, data_source='EBI-GOA non-IntAct')
                    session.add(source)
        session.commit()
        print(session.query(Interaction).count())

