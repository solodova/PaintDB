import csv
from Schema1 import Interactor, Interaction, OrthologEcoli, InteractionReference, InteractionSource
from Parsers.Parser import is_experimental_psimi

def parse_ecoli_dip(session):
    with open('Ecoli/DIP.txt') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
            interactors = []

            ids_A = row['ID interactor A'].split('|')
            ids_B = row['ID interactor B'].split('|')
            refseq_A, uniprotkb_A, refseq_B, uniprotkb_B = '', '', '', ''
            for id in ids_A:
                fields = id.split(':')
                if fields[0] == 'refseq':
                    refseq_A = fields[1]
                elif fields[0] == 'uniprotkb':
                    uniprotkb_A = fields[1]
            for id in ids_B:
                fields = id.split(':')
                if fields[0] == 'refseq':
                    refseq_B = fields[1]
                elif fields[0] == 'uniprotkb':
                    uniprotkb_B = fields[1]

            orthologs_A, orthologs_B = [], []
            if uniprotkb_A != '':
                orthologs_A = session.query(OrthologEcoli).filter(OrthologEcoli.ortholog_uniprot == uniprotkb_A).all()
            if (len(orthologs_A) == 0) & (refseq_A != ''):
                orthologs_A = session.query(OrthologEcoli).filter(OrthologEcoli.ortholog_refseq == refseq_A).all()
            if uniprotkb_B != '':
                orthologs_B = session.query(OrthologEcoli).filter(OrthologEcoli.ortholog_uniprot == uniprotkb_B).all()
            if (len(orthologs_B) == 0) & (refseq_B != ''):
                orthologs_B = session.query(OrthologEcoli).filter(OrthologEcoli.ortholog_refseq == refseq_B).all()

            for ortholog_A in orthologs_A:
                for ortholog_B in orthologs_B:
                    if (ortholog_A is not None) and (ortholog_B is not None):
                        if ortholog_A.strain_protein == ortholog_B.strain_protein:
                            interactors.append([[ortholog_A.protein, ortholog_A.ortholog_id], [ortholog_B.protein,
                                                                                               ortholog_B.ortholog_id]])

            for interactor_pair in interactors:
                is_new = 0
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
                    is_new = 1
                    interaction = Interaction(strain=interactor_pair[0][0].strain, interactors=[interactor_pair[0][0],
                                              interactor_pair[1][0]], type='p-p', ortholog_derived='fe')
                    session.add(interaction), session.commit()

                detections, pmids, types, list = [], [], [], []
                if row['Interaction detection method(s)'] != '-':
                    detections = row['Interaction detection method(s)'].split('|')
                    list.append(detections)
                if row['Publication Identifier(s)'] != '-':
                    pmids = row['Publication Identifier(s)'].split('|')
                    list.append(pmids)
                if row['Interaction type(s)'] != '-':
                    types = row['Interaction type(s)'].split('|')
                    list.append(types)

                interactor_a, interactor_b = '', ''
                if interaction.interactors[0] == interactor_pair[0][0]:
                    interactor_a = interactor_pair[0][1]
                    interactor_b = interactor_pair[1][1]
                else:
                    interactor_b = interactor_pair[0][1]
                    interactor_a = interactor_pair[1][1]

                for num in range(0, len(list[0])):
                    type = types[num].split('(')[1][:-1]
                    pmid = pmids[num * 2].split('pubmed:')[1]
                    detection = detections[num].split('(')[1][:-1]
                    # there are more than one pmid sometimes
                    reference = InteractionReference(interaction_id=interaction.id,
                                                     detection_method=detection,
                                                     pmid=pmid,
                                                     source_db=row['Source database(s)'].split('(')[1][:-1],
                                                     interactor_a = interactor_a,
                                                     interactor_b = interactor_b)
                    session.add(reference)

                    if is_new:
                        if interaction.is_experimental is None:
                            if is_experimental_psimi(row['Interaction detection method(s)'].split('MI:')[1][:4]):
                                interaction.is_experimental = 1
                            else:
                                interaction.is_experimental = 0
                        elif is_experimental_psimi(row['Interaction detection method(s)'].split('MI:')[1][:4]):
                            interaction.is_experimental = 1

                source = session.query(InteractionSource).filter(InteractionSource.interaction_id == interaction.id,
                                                                 InteractionSource.data_source == 'DIP').first()

                if source is None:
                    source = InteractionSource(interaction_id=interaction.id, data_source='DIP')
                    session.add(source)
        session.commit()
        print(session.query(Interaction).count())