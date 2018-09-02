import csv
from Schema1 import Interactor, OrthologPseudomonas, Interaction, InteractionReference, InteractionSource

inparalogs = {}

def parse(session):
    get_inparalogs('Data/Ortholog/PAO1-PA14.csv')
    parse_orthologs('Data/Ortholog/PAO1-PA14.csv', session)
    parse_ortholog_interactions(session)

def get_inparalogs(ortholog_file):
    with open(ortholog_file) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if row['Strain 1 Inparalogs (Locus Tag/Name)'] != '':
                strain1_inparalogs = row['Strain 1 Inparalogs (Locus Tag/Name)'].split(']')

                for inparalog in strain1_inparalogs:
                    if inparalog == '': continue
                    inparalog_id=inparalog.split('[')[0]
                    if inparalog_id[0] == ';':
                        inparalog_id = inparalog_id[1:]
                    if inparalog_id in inparalogs:
                        inparalogs[inparalog_id].append(row['Locus Tag (Strain 2)'])
                    else:
                        inparalogs[inparalog_id] = [row['Locus Tag (Strain 2)']]

            if row['Strain 2 Inparalogs (Locus Tag/Name)'] != '':
                strain2_inparalogs = row['Strain 2 Inparalogs (Locus Tag/Name)'].split(']')

                for inparalog in strain2_inparalogs:
                    if inparalog == '': continue
                    inparalog_id = inparalog.split('[')[0]
                    if inparalog_id[0] == ';':
                        inparalog_id = inparalog_id[1:]
                    if row['Locus Tag (Strain 1)'] in inparalogs:
                        inparalogs[row['Locus Tag (Strain 1)']].append(inparalog_id)
                    else:
                        inparalogs[row['Locus Tag (Strain 1)']] = [inparalog_id]

def parse_orthologs(ortholog_file, session):
    with open(ortholog_file) as csvfile:
        reader = csv.DictReader(csvfile)
        orthologs = []
        for row in reader:
            if row['Locus Tag (Strain 1)'] in inparalogs:
                if row['Locus Tag (Strain 2)'] in inparalogs[row['Locus Tag (Strain 1)']]:
                    continue

            if session.query(Interactor).get(row['Locus Tag (Strain 1)']) is not None:
                if session.query(Interactor).get(row['Locus Tag (Strain 2)']) is not None:

                    ortholog1 = OrthologPseudomonas(protein_id=row['Locus Tag (Strain 1)'], strain_protein='PAO1',
                                                    ortholog_id=row['Locus Tag (Strain 2)'], strain_ortholog='PA14')
                    ortholog2 = OrthologPseudomonas(protein_id=row['Locus Tag (Strain 2)'], strain_protein='PA14',
                                                    ortholog_id=row['Locus Tag (Strain 1)'], strain_ortholog='PAO1')
                    orthologs.append(ortholog1), orthologs.append(ortholog2)
        session.add_all(orthologs)
        session.commit()


def parse_ortholog_interactions(session):
    all_interactions =session.query(Interaction).all()
    for interaction in all_interactions:
        interactor_pairs, ortholog_interactors = [], [[], []]
        num = 0
        for interactor in interaction.interactors:
            if interactor.type == 'p':
                for ortholog in interactor.pseudomonas_orthologs:
                    if ortholog is not None:
                        ortholog_interactors[num].append(ortholog.protein)
            else:
                ortholog_interactors[num].append(interactor)
            num += 1

        for interactor1 in ortholog_interactors[0]:
            for interactor2 in ortholog_interactors[1]:
                interactor_pairs.append([interactor1, interactor2])

        for interactor_pair in interactor_pairs:

            homogenous = (interactor_pair[0] == interactor_pair[1])
            new_interaction = session.query(Interaction).filter(Interaction.interactors.contains(interactor_pair[0]),
                                                                Interaction.interactors.contains(interactor_pair[1]),
                                                                Interaction.homogenous == homogenous).first()
            # should you add new references if the interaction already exists?
            if new_interaction is None:
                strain = 'PAO1'
                if interaction.strain == 'PAO1':
                    strain = 'PA14'
                new_interaction = Interaction(strain=strain, type=interaction.type, interactors=interactor_pair,
                                              homogenous=homogenous)
                session.add(new_interaction), session.commit()

            for reference in interaction.references:
                new_ref = session.query(InteractionReference).filter_by(psimi_detection = reference.psimi_detection,
                                                                        detection_method=reference.detection_method,
                                                                        author_ln=reference.author_ln,
                                                                        pub_date=reference.pub_date,
                                                                        pmid=reference.pmid,
                                                                        psimi_type=reference.psimi_type,
                                                                        interaction_type=reference.interaction_type,
                                                                        psimi_db= reference.psimi_db,
                                                                        source_db=reference.source_db,
                                                                        confidence=reference.confidence,
                                                                        comment = reference.comment,
                                                                        interactor_a=interaction.interactors[0].id,
                                                                        interactor_b=interaction.interactors[1].id).first()

                if new_ref is None:
                    new_ref = InteractionReference(psimi_detection=reference.psimi_detection,
                                                         detection_method=reference.detection_method,
                                                         author_ln=reference.author_ln, pub_date=reference.pub_date,
                                                         pmid=reference.pmid, psimi_type = reference.psimi_type,
                                                         interaction_type=reference.interaction_type,
                                                         psimi_db=reference.psimi_db, source_db=reference.source_db,
                                                         confidence=reference.confidence, comment = reference.comment,
                                                         interactor_a = interaction.interactors[0].id,
                                                         interactor_b = interaction.interactors[1].id)
                    new_interaction.references.append(new_ref)
                elif new_ref not in new_interaction.references:
                    new_interaction.references.append(new_ref)

            for source in interaction.sources:
                new_source = session.query(InteractionSource).filter_by(
                    data_source = source.data_source, is_experimental = source.is_experimental).first()
                if new_source is None:
                    new_source = InteractionSource(data_source = source.data_source,
                                                    is_experimental = source.is_experimental)
                    new_interaction.sources.append(new_source)
                elif new_source not in new_interaction.sources:
                    new_interaction.sources.append(new_source)

        session.commit()
