import csv
from Schema1 import Interactor, Metabolite, Interaction, InteractionReference, OrthologEcoli, InteractionXref, InteractorXref, \
    InteractionSource
from Schema1 import is_experimental_psimi

cols = ['interactor_A', 'interactor_B', 'altID_A', 'altID_B', 'alias_A', 'alias_B', 'detection', 'publication',
        'publication_ID', 'taxid_A', 'taxid_B', 'type', 'source_db', 'identifier', 'confidence']

def parse_ecoli_bindingdb(session):
    with open('Data/Ecoli/PSICQUIC/BindingDB.txt') as csvfile:
        reader = csv.DictReader(csvfile, fieldnames=cols, delimiter = '\t')

        # iterate through each interaction
        for row in reader:
            #if (row['interactor_A'] == '-') | (row['interactor_B'] == '-'): continue
            uniprot_A, refseq_A, orthologs_A = None, None, None
            uniprot_B, refseq_B, orthologs_B = None, None, None
            pubchem, chebi = None, None
            metabolite_info, metabolite, orthologs = None, None, None

            if 'uniprotkb' in row['interactor_A']:
                uniprot_A = row['interactor_A'].split('uniprotkb:')[1].split('|')[0]
            if 'refseq' in row['interactor_A']:
                refseq_A = row['interactor_A'].split('refseq:')[1].split('|')[0]

            if uniprot_A is not None:
                orthologs_A = session.query(OrthologEcoli).filter(OrthologEcoli.ortholog_uniprot == uniprot_A).all()
            if (orthologs_A is None) and (refseq_A is not None):
                orthologs_A = session.query(OrthologEcoli).filter(OrthologEcoli.ortholog_refseq == refseq_A).all()
            if orthologs_A is None & ((uniprot_A is not None) | (refseq_A is not None)): continue

            if 'uniprotkb' in row['interactor_B']:
                uniprot_B = row['interactor_B'].split('uniprotkb:')[1].split('|')[0]
            if 'refseq' in row['interactor_B']:
                refseq_B = row['interactor_B'].split('refseq:')[1].split('|')[0]

            if uniprot_B is not None:
                orthologs_B = session.query(OrthologEcoli).filter(OrthologEcoli.ortholog_uniprot == uniprot_B).all()
            if (orthologs_B is None) and (refseq_B is not None):
                orthologs_B = session.query(OrthologEcoli).filter(OrthologEcoli.ortholog_refseq == refseq_B).all()
            if orthologs_B is None & ((uniprot_B is not None) | (refseq_B is not None)): continue

            if (orthologs_A is None) and (orthologs_B is None): continue

            if orthologs_A is None:
                if 'chebi' in row['interactor_A']:
                    chebi = row['interactor_A'].split('CHEBI:')[1].split('|')[0][:-1]
                if 'pubchem' in row['altID_A']:
                    pubchem = row['altID_A'].split('pubchem:')[1].split('|')[0]
                if (chebi is None) & ('chebi' in row['altID_A']):
                    chebi = row['altID_A'].split('CHEBI:')[1].split('|')[0][:-1]
                if (chebi is None) & (pubchem is None): continue
                orthologs = orthologs_B
            elif orthologs_B is None:
                if 'chebi' in row['interactor_B']:
                    chebi = row['interactor_B'].split('CHEBI:')[1].split('|')[0][:-1]
                if 'pubchem' in row['altID_B']:
                    pubchem = row['altID_B'].split('pubchem:')[1].split('|')[0]
                if (chebi is None) & ('chebi' in row['altID_B']):
                    chebi = row['altID_B'].split('CHEBI:')[1].split('|')[0][:-1]
                if (chebi is None) & (pubchem is None): continue
                orthologs = orthologs_A

            if (chebi is not None) | (pubchem is not None):
                id = None
                if chebi is not None:
                    id = chebi
                    metabolite = session.query(Metabolite).filter(Metabolite.chebi == chebi).first()
                if (metabolite is None) & (pubchem is not None):
                    id = pubchem
                    metabolite = session.query(Metabolite).filter(Metabolite.pubchem == pubchem).first()
                if metabolite is None:
                    metabolite = Metabolite(id = id, chebi=chebi, pubchem=pubchem)
                    session.add(metabolite), session.commit()

            interactors = []
            if metabolite is None:
                for ortholog_A in orthologs_A:
                    for ortholog_B in orthologs_B:
                        if (ortholog_A is not None) and (ortholog_B is not None):
                            if ortholog_A.strain_protein == ortholog_B.strain_protein:
                                interactors.append([[ortholog_A.protein, ortholog_A.ortholog_id],
                                                    [ortholog_B.protein, ortholog_B.ortholog_id]])
            else:
                for ortholog in orthologs:
                    interactors.append([[metabolite, metabolite.id], [ortholog.protein, ortholog.ortholog_id]])

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
                else:
                    interaction = Interaction(strain=interactor_pair[0][0].strain,
                                              interactors=[interactor_pair[0][0], interactor_pair[1][0]],
                                              type=(interactor_pair[0][0].type + '-' + interactor_pair[1][0]),
                                              ortholog_derived='fe')
                    session.add(interaction), session.commit()

                interactor_a, interactor_b = None, None
                if interaction.interactors[0] == interactor_pair[0][0]:
                    interactor_a = interactor_pair[0][1]
                    interactor_b = interactor_pair[1][1]
                else:
                    interactor_b = interactor_pair[0][1]
                    interactor_a = interactor_pair[1][1]

                psimi_detection, psimi_db, psimi_type, author, date, confidences = None, None, None, None, None, [None]
                    if 'MI' in row['Interaction detection method(s)']:
                        psimi_detection = row['Interaction detection method(s)'].split('MI:')[1][:4]
                    if 'MI' in row['Interaction type(s)']:
                        psimi_type = row['Interaction type(s)'].split('MI:')[1][:4]
                    if 'MI' in row['Source database(s)']:
                        psimi_db = row['Source database(s)'].split('MI:')[1][:4]
                    if row['Publication 1st author(s)'] != '-':
                        author = row['Publication 1st author(s)'].split(' ')[0]
                        date = row['Publication 1st author(s)'].split('(')[1][:-1]
                    if ('intact-miscore' in row['Confidence value(s)']) | (
                            'author score' in row['Confidence value(s)']):
                        del confidences[0]
                        confidence_ids = row['Confidence value(s)'].split('|')
                        for confidence in confidence_ids:
                            if (confidence.split(':')[0] == 'intact-miscore') | \
                                    (confidence.split(':')[0] == 'author score'):
                                confidences.append(confidence)
                    for confidence in confidences:
                        reference = InteractionReference(interaction_id=interaction.id,
                                                         psimi_detection=psimi_detection,
                                                         detection_method=
                                                         row['Interaction detection method(s)'].split('(')[1][:-1],
                                                         author_ln=author,
                                                         pub_date=date,
                                                         pmid=
                                                         row['Publication Identifier(s)'].split('pubmed:')[1].split(
                                                             '|')[0],
                                                         psimi_type=psimi_type,
                                                         interaction_type=row['Interaction type(s)'].split('(')[1][:-1],
                                                         psimi_db=psimi_db,
                                                         source_db=row['Source database(s)'].split('(')[1][:-1],
                                                         confidence=confidence,
                                                         interactor_a_id=interactor_a,
                                                         interactor_b_id=interactor_b)
                        session.add(reference)

                    source = session.query(InteractionSource).filter(
                        InteractionSource.interaction_id == interaction.id,
                        InteractionSource.data_source == 'IMEx').first()

                    if source is None:
                        source = InteractionSource(interaction_id=interaction.id, data_source='IMEx')
                        session.add(source)
            session.commit()
            print(session.query(Interaction).count())