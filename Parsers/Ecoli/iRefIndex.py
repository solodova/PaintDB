import csv
from Schema1 import OrthologEcoli, Interactor, Interaction, InteractionReference

def parse_Ecoli_iRefIndex(session):
    with open('Ecoli/PSICQUIC/iRefIndex.txt') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
            if (row['#ID(s) interactor A'] == '-') | (row['ID(s) interactor B'] == '-'): continue
            interactors, orthologs_A, orthologs_B = [], [], []
            if row['#ID(s) interactor A'].split(':')[0] == 'uniprotkb':
                orthologs_A = session.query(OrthologEcoli).filter(
                    OrthologEcoli.ortholog_uniprot == row['#ID(s) interactor A'].split(':')[1]).all()
            elif row['#ID(s) interactor A'].split(':')[0] == 'refseq':
                orthologs_A = session.query(OrthologEcoli).filter(
                    OrthologEcoli.ortholog_refseq == row['#ID(s) interactor A'].split(':')[1]).all()
            if row['ID(s) interactor B'].split(':')[0] == 'uniprotkb':
                orthologs_B = session.query(OrthologEcoli).filter(
                    OrthologEcoli.ortholog_uniprot == row['ID(s) interactor B'].split(':')[1]).all()
            elif row['ID(s) interactor B'].split(':')[0] == 'refseq':
                orthologs_B = session.query(OrthologEcoli).filter(
                    OrthologEcoli.ortholog_refseq == row['ID(s) interactor B'].split(':')[1]).all()

            for ortholog_A in orthologs_A:
                if ortholog_A != None:
                    interactor_A = session.query(Interactor).filter(Interactor.id == ortholog_A.protein_id).one()
                    for ortholog_B in orthologs_B:
                        if ortholog_B != None:
                            interactor_B = session.query(Interactor).filter(
                                Interactor.id == ortholog_B.protein_id).one()
                            if (interactor_A.strain == interactor_B.strain):
                                interactors.append([interactor_A, interactor_B])

            for interactor_pair in interactors:
                homogenous = (interactor_pair[0] == interactor_pair[1])
                interaction = session.query(Interaction).filter(
                    (Interaction.interactors.contains(interactor_pair[0])),
                    (Interaction.interactors.contains(interactor_pair[1])),
                    (Interaction.homogenous == homogenous)).first()
                if (interaction != None):
                    if (interaction.ortholog_derived == None):
                        interaction.ortholog_derived = 'confirmed from E.coli'
                    elif ('from E. coli' not in interaction.ortholog_derived):
                        interaction.ortholog_derived += ', confirmed from E. coli'
                    session.commit()
                if (interaction == None):
                    interaction = Interaction(strain=interactor_pair[0].strain, interactors=interactor_pair,
                                              type=(interactor_pair[0].type + '-' + interactor_pair[1].type),
                                              is_experimental=0, ortholog_derived='from E. coli')
                    session.add(interaction), session.commit()
                author, date, ref, type, pmid, detection = None, None, None, None, None, None
                if (row['Publication 1st author(s)'] != '-'):
                    author = row['Publication 1st author(s)'].split('-')[0]
                    date = row['Publication 1st author(s)'].split('-')[1]
                    ref = row['Publication 1st author(s)']
                if (row['Interaction type(s)'] != '-'):
                    type = row['Interaction type(s)'].split('(')[1][:-1]
                if (row['Publication Identifier(s)'] != '-'):
                    pmid = row['Publication Identifier(s)'].split('med:')[1][:8]
                if (row['Interaction detection method(s)'] != '-'):
                    detection = row['Interaction detection method(s)'].split('(')[1][:-1]
                # there are more than one pmid sometimes
                reference = InteractionReference(interaction_type=type, interaction_id=interaction.id,
                                                 confidence_score=row['Confidence value(s)'], pmid=pmid,
                                                 detection_method=detection, author_last_name=author,
                                                 publication_date=date, publication_ref=ref,
                                                 source_db=row['Source database(s)'].split('(')[1][:-1],
                                                 interactor_a_id=row['#ID(s) interactor A'].split(':')[1],
                                                 interactor_b_id=row['ID(s) interactor B'].split(':')[1])
                session.add(reference)
        session.commit()
        print(session.query(Interaction).count())