import csv
from Schema1 import Interactor, OrthologPseudomonas, Interaction, InteractionReference, InteractionSource
from sqlalchemy import or_

inparalogs = {}

def parse_ortholuge(session):
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
                    if inparalog_id in inparalogs:
                        inparalogs[inparalog_id].append(row['Locus Tag (Strain 1)'])
                    else:
                        inparalogs[inparalog_id] = [row['Locus Tag (Strain 1)']]

def parse_orthologs(ortholog_file, session):
    with open(ortholog_file) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if row['Locus Tag (Strain 1)'] in inparalogs:
                if row['Locus Tag (Strain 2)'] in inparalogs[row['Locus Tag (Strain 1)']]:
                    continue
            if row['Locus Tag (Strain 2)'] in inparalogs:
                if row['Locus Tag (Strain 1)'] in inparalogs[row['Locus Tag (Strain 2)']]:
                    continue

            if session.query(Interactor).filter(Interactor.id == row['Locus Tag (Strain 1)']).first() is not None:
                if session.query(Interactor).filter(Interactor.id == row['Locus Tag (Strain 2)']).first() is not None:
                    ortholog1 = OrthologPseudomonas(protein_id=row['Locus Tag (Strain 1)'], strain_protein='PAO1',
                                                    ortholog_id=row['Locus Tag (Strain 2)'], strain_ortholog='PA14')
                    ortholog2 = OrthologPseudomonas(protein_id=row['Locus Tag (Strain 2)'], strain_protein='PA14',
                                                    ortholog_id=row['Locus Tag (Strain 1)'], strain_ortholog='PAO1')
                    session.add(ortholog1), session.add(ortholog2)
        session.commit()


def parse_ortholog_interactions(session):
    all_interactions =session.query(Interaction).filter(or_(Interaction.type == 'p-p',
                                                         Interaction.type == 'm-p',
                                                         Interaction.type == 'p-m')).all()
    for interaction in all_interactions:
        new_interactors, ortholog_interactors = [], [[], []]
        num = 0
        ortho = []
        for interactor in interaction.interactors:
            if interactor.type == 'p':
                for ortholog in interactor.pseudomonas_orthologs:
                    if ortholog is not None:
                        ortholog_interactors[num].append(session.query(Interactor).filter(
                            Interactor.id == ortholog.ortholog_id).one())
            else:
                ortholog_interactors[num].append(interactor)
            num += 1

        for interactor1 in ortholog_interactors[0]:
            for interactor2 in ortholog_interactors[1]:
                new_interactors.append([interactor1, interactor2])

        ortho.append(new_interactors)

        for interactor_pair in new_interactors:

            homogenous = (interactor_pair[0] == interactor_pair[1])
            new_interaction = session.query(Interaction).filter(Interaction.interactors.contains(interactor_pair[0]),
                                                                Interaction.interactors.contains(interactor_pair[1]),
                                                                Interaction.homogenous == homogenous).first()
            # should you add new references if the interaction already exists?
            if new_interaction is not None:
                if new_interaction.ortholog_derived is None:
                    new_interaction.ortholog_derived = 'cf' + interaction.strain
                elif ('f' + interaction.strain) not in new_interaction.ortholog_derived:
                    new_interaction.ortholog_derived += ', cf' + interaction.strain
                session.commit()
            else:
                strain = 'PAO1'
                if interaction.strain == 'PAO1':
                    strain = 'PA14'
                new_interaction = Interaction(strain=strain, type=interaction.type,
                                              ortholog_derived='f' + interaction.strain, interactors=interactor_pair,
                                              homogenous=homogenous)
                session.add(new_interaction), session.commit()

            for reference in interaction.references:
                new_ref = session.query(InteractionReference).filter(
                    InteractionReference.psimi_detection == reference.psimi_detection,
                    InteractionReference.detection_method==reference.detection_method,
                    InteractionReference.author_ln==reference.author_ln,
                    InteractionReference.pub_date==reference.pub_date,
                    InteractionReference.pmid==reference.pmid,
                    InteractionReference.psimi_type==reference.psimi_type,
                    InteractionReference.interaction_type==reference.interaction_type,
                    InteractionReference.psimi_db== reference.psimi_db,
                    InteractionReference.source_db==reference.source_db,
                    InteractionReference.confidence==reference.confidence,
                    InteractionReference.comment == reference.comment,
                    InteractionReference.interactor_a==interaction.interactors[0].id,
                    InteractionReference.interactor_b==interaction.interactors[1].id).first()

                if new_ref is None:
                    new_reference = InteractionReference(psimi_detection=reference.psimi_detection,
                                                         detection_method=reference.detection_method,
                                                         author_ln=reference.author_ln,
                                                         pub_date=reference.pub_date,
                                                         pmid=reference.pmid,
                                                         psimi_type = reference.psimi_type,
                                                         interaction_type=reference.interaction_type,
                                                         psimi_db=reference.psimi_db,
                                                         source_db=reference.source_db,
                                                         confidence=reference.confidence,
                                                         comment = reference.comment,
                                                         interactor_a = interaction.interactors[0].id,
                                                         interactor_b = interaction.interactors[1].id)
                    session.add(new_reference), session.commit()
                    new_interaction.references.append(new_reference)

            for source in interaction.sources:
                new_source = session.query(InteractionSource).filter(
                    InteractionSource.data_source == source.data_source,
                    InteractionSource.is_experimental == source.is_experimental).first()
                if new_source is None:
                    new_source = InteractionSource(data_source = source.data_source,
                                                    is_experimental = source.is_experimental)
                    session.add(new_source), session.commit()
                    new_interaction.sources.append(new_source)
                elif new_source not in new_interaction.sources:
                    new_interaction.sources.append(new_source)

        session.commit()
