import csv
from Schema import Interactor, OrthologPseudomonas, Interaction, InteractionReference, InteractionSource

# dicts of form {strain 1 id: [strain 2 inparalogs]}
inparalogs = {}

def parse(session):
    parse_inparalogs('Data/Ortholog/PAO1-PA14.csv')
    parse_orthologs('Data/Ortholog/PAO1-PA14.csv', session)
    parse_ortholog_interactions(session)

def parse_inparalogs(file):
    # get all inparalogs and store in inparalogs dict
    with open(file) as csvfile:
        reader = csv.DictReader(csvfile)

        for row in reader:
            # ignore Non/Borderline SSD inparalogs (they will be ignored when parsing orthologs anyway)
            if (row['Ortholuge Class'] != '') & (row['Ortholuge Class'] != 'SSD'): continue

            #check for strain 1 (PAO1) inparalogs
            if row['Strain 1 Inparalogs (Locus Tag/Name)'] != '':
                strain1_inparalogs = row['Strain 1 Inparalogs (Locus Tag/Name)'].split(']')

                # for all inparalogs, add them to inparalog dict
                for inparalog in strain1_inparalogs:
                    if inparalog == '': continue
                    inparalog_id = inparalog.split('[')[0]
                    if inparalog_id[0] == ';':
                        inparalog_id = inparalog_id[1:]
                    if inparalog_id in dict:
                        inparalogs[inparalog_id].append(row['Locus Tag (Strain 2)'])
                    else:
                        inparalogs[inparalog_id] = [row['Locus Tag (Strain 2)']]

            # check for strain 2 (PA14) inparalogs
            if row['Strain 2 Inparalogs (Locus Tag/Name)'] != '':
                strain2_inparalogs = row['Strain 2 Inparalogs (Locus Tag/Name)'].split(']')

                for inparalog in strain2_inparalogs:
                    if inparalog == '': continue
                    inparalog_id = inparalog.split('[')[0]
                    if inparalog_id[0] == ';':
                        inparalog_id = inparalog_id[1:]
                    if row['Locus Tag (Strain 1)'] in dict:
                        inparalogs[row['Locus Tag (Strain 1)']].append(inparalog_id)
                    else:
                        inparalogs[row['Locus Tag (Strain 1)']] = [inparalog_id]

def parse_orthologs(ortholog_file, session):
    with open(ortholog_file) as csvfile:
        reader = csv.DictReader(csvfile)
        orthologs = []
        for row in reader:
            # ignore inparalogs
            if row['Locus Tag (Strain 1)'] in inparalogs:
                if row['Locus Tag (Strain 2)'] in inparalogs[row['Locus Tag (Strain 1)']]:
                    continue
            if session.query(Interactor).get(row['Locus Tag (Strain 1)']) is not None:
                if session.query(Interactor).get(row['Locus Tag (Strain 2)']) is not None:
                    # add a PA14 ortholog to the PAO1 ortholog partner, and vice versa
                    ortholog1 = OrthologPseudomonas(protein_id=row['Locus Tag (Strain 1)'], strain_protein='PAO1',
                                                    ortholog_id=row['Locus Tag (Strain 2)'], strain_ortholog='PA14')
                    ortholog2 = OrthologPseudomonas(protein_id=row['Locus Tag (Strain 2)'], strain_protein='PA14',
                                                    ortholog_id=row['Locus Tag (Strain 1)'], strain_ortholog='PAO1')
                    orthologs.append(ortholog1), orthologs.append(ortholog2)
        session.add_all(orthologs)
        session.commit()


#function to parse all ortholog interactions
def parse_ortholog_interactions(session):
    # query for all current interactions from PAO1 and PA14 sources
    all_interactions = session.query(Interaction).all()
    # iterate through each interaction, see if interactors have orthologs, and create new interactions in
    # other strain if they do
    for interaction in all_interactions:
        # ortholog interactors is interactors from opposite strain from that in interaction
        interactor_pairs, ortholog_interactors = [], [[], []]
        num = 0
        for interactor in interaction.interactors:
            # if the interactor is a protein, add its pseudomonas orthologs to ortholog_interactors[num]
            if interactor.type == 'p':
                for ortholog in interactor.pseudomonas_orthologs:
                    if ortholog is not None:
                        # add the interactor's psuedomonas ortholog to ortholog_interactors
                        # also add the interactor id for creation of interaction reference later
                        ortholog_interactor = session.query(Interactor).get(ortholog.ortholog_id)
                        ortholog_interactors[num].append([ortholog_interactor, interactor.id])
            # if the interactor is a metabolite, add it as is
            else:
                ortholog_interactors[num].append([interactor, interactor.id])
            num += 1

        # create interactor pairs from ortholog interactors
        for interactor1 in ortholog_interactors[0]:
            for interactor2 in ortholog_interactors[1]:
                interactor_pairs.append([interactor1, interactor2])

        # iterate through each interactor pair, create interaction if it doesnt already exist
        for interactor_pair in interactor_pairs:
            homogenous = (interactor_pair[0] == interactor_pair[1])
            new_interaction = session.query(Interaction).filter(Interaction.interactors.contains(interactor_pair[0]),
                                                                Interaction.interactors.contains(interactor_pair[1]),
                                                                Interaction.homogenous == homogenous).first()
            if new_interaction is None:
                # set strain for new interaction to opposite of original interaction
                strain = 'PAO1'
                if interaction.strain == 'PAO1':
                    strain = 'PA14'
                # set ortholog derived to the original interaction strain
                new_interaction = Interaction(strain=strain, type=interaction.type, interactors=interactor_pair,
                                              homogenous=homogenous, ortholog_derived = interaction.strain)
                session.add(new_interaction), session.commit()

            # in case the interaction already existed, make sure interactor_a and interactor_b variables for
            # new interaction reference match up with the first and second interactors of the existing
            # interaction (so it's easy to see which Pseudomonas interactor matches up with which Ecoli
            # ortholog)
            interactor_a, interactor_b = None, None
            if new_interaction.interactors[0] == interactor_pair[0]:
                interactor_a = interactor_pair[0][1]
                interactor_b = interactor_pair[1][1]
            else:
                interactor_b = interactor_pair[0][1]
                interactor_a = interactor_pair[1][1]

            # iterate through all of the original interaction references, create new reference
            # with same fields except with added interactor_a and interactor_b attributes to show original
            # interactors from which interaction was derived from
            for reference in interaction.references:
                new_ref = session.query(InteractionReference).filter_by(detection_method=reference.detection_method,
                                                                        author_ln=reference.author_ln,
                                                                        pub_date=reference.pub_date,
                                                                        pmid=reference.pmid,
                                                                        interaction_type=reference.interaction_type,
                                                                        source_db=reference.source_db,
                                                                        confidence=reference.confidence,
                                                                        comment = reference.comment,
                                                                        interactor_a=interactor_a,
                                                                        interactor_b=interactor_b).first()

                if new_ref is None:
                    # if the new_ref doesn't exist, create and add it to the new interaction's reference list
                    # and add the original reference's sources to the new ones sources
                    new_ref = InteractionReference(detection_method=reference.detection_method,
                                                   author_ln=reference.author_ln, pub_date=reference.pub_date,
                                                   pmid=reference.pmid, interaction_type=reference.interaction_type,
                                                   source_db=reference.source_db,confidence=reference.confidence,
                                                   comment = reference.comment,
                                                   interactor_a = interactor_a, interactor_b = interactor_b)
                    new_interaction.references.append(new_ref)
                    new_ref.sources = reference.sources
                else:
                    # if the new reference did exist, add the new interaction and original interactions sources
                    # to new reference's attributes if they were not there already
                    if new_interaction not in new_ref.interactions:
                        new_ref.interactions.append(new_interaction)
                    for source in reference.sources:
                        if source is not None:
                            if source not in new_ref.sources:
                                new_ref.sources.append(source)

            # for each source in the original interaction's sources, add it to the new interaction's source list if
            # it isn't already there
            for source in interaction.sources:
                if source not in new_interaction.sources:
                    new_interaction.sources.append(source)

    session.commit()
    print('p_orthologs', session.query(Interaction).count())