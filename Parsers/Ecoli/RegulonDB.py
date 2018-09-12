import csv
from Schema import OrthologEcoli, Interaction, InteractionReference, InteractionSource

def parse(session):
    with open('Data/Ecoli/RegulonDB.csv') as csvfile:
        reader = csv.DictReader(csvfile)
        # since all the interactions from here will use the same source, create and add it at the beginning
        # Note: since no references are available, is_experimental is set to 2
        source = InteractionSource(data_source = 'RegulonDB(Ecoli)', is_experimental = 2)
        session.add(source), session.commit()

        for row in reader:
            interactors = []

            orthologs_A = session.query(OrthologEcoli).filter_by(
                ortholog_name = (row['TF name'][0].lower() + row['TF name'][1:])).all()
            # if no orthologs for first interactor were found, skip to next interaction
            if orthologs_A is None: continue
            orthologs_B = session.query(OrthologEcoli).filter_by(ortholog_name = row['Regulated gene']).all()
            # if no orthologs for second interactor were found, skip to next interaction
            if orthologs_B is None: continue

            # iterate through each ortholog in ortholog A and B to create interactor pairs from their
            # respective pseudomonas proteins
            for ortholog_A in orthologs_A:
                for ortholog_B in orthologs_B:
                    # only add the pseudomonas interactors if their strains match
                    if ortholog_A.strain_protein == ortholog_B.strain_protein:
                        # make sure to add ortholog id for creating the interaction reference later
                        interactors.append([[ortholog_A.protein, ortholog_A.ortholog_id],
                                            [ortholog_B.protein, ortholog_B.ortholog_id]])

            # iterate through each interactor pair, create a new interaction if it doesnt exist yet
            for interactor_pair in interactors:
                homogenous = (interactor_pair[0][0] == interactor_pair[1][0])
                interaction = session.query(Interaction).filter(Interaction.interactors.contains(interactor_pair[0][0]),
                                                                Interaction.interactors.contains(interactor_pair[1][0]),
                                                                Interaction.homogenous == homogenous).first()

                if interaction is None:
                    # if interaction is None, make ortholog_derived = Ecoli and add source to interaction sources
                    interaction = Interaction(strain=interactor_pair[0][0].strain,
                                              interactors=[interactor_pair[0][0], interactor_pair[1][0]], type='p-p',
                                              ortholog_derived='Ecoli')
                    interaction.sources.append(source)
                    session.add(interaction), session.commit()
                elif source not in interaction.sources:
                    interaction.sources.append(source)

                # in case the interaction already existed, make sure interactor_a and interactor_b variables for
                # new interaction reference match up with the first and second interactors of the existing
                # interaction (so it's easy to see which Pseudomonas interactor matches up with which Ecoli
                # ortholog)
                interactor_a, interactor_b = None, None
                if interaction.interactors[0] == interactor_pair[0][0]:
                    interactor_a = interactor_pair[0][1]
                    interactor_b = interactor_pair[1][1]
                else:
                    interactor_b = interactor_pair[0][1]
                    interactor_a = interactor_pair[1][1]

                type= 'TF/sigma-binding site (' + row['Regulatory effect'] + 'regulation)'
                comment= interactor_pair[0][1] + ' regulates(' +row['Regulatory effect'] + ') ' + interactor_pair[1][1]
                # create a reference for each evidence type listed for interaction
                for evidence in row['Evidence'][1:-1].split(', '):
                    # check if interaction reference already exists in db
                    reference = session.query(InteractionReference).filter_by(detection_method = evidence,
                                                                              interaction_type = type,
                                                                              source_db = 'regulondb',
                                                                              confidence = row['Evidence type'],
                                                                              comment = comment,
                                                                              interactor_a = interactor_a,
                                                                              interactor_b = interactor_b).first()
                    if reference is None:
                        # if reference is None, add reference to interaction references list and add source
                        # to reference sources list
                        reference = InteractionReference(detection_method = evidence, interaction_type = type,
                                                         comment=comment, source_db='regulondb',
                                                         confidence = row['Evidence type'],
                                                         interactor_a=interactor_a, interactor_b=interactor_b)
                        interaction.references.append(reference)
                        reference.sources.append(source)
                    # if reference exists, check that its interactions contains interaction, and sources contains
                    # source, and add if they are not present
                    else:
                        if interaction not in reference.interactions:
                            interaction.references.append(reference)
                        if source not in reference.sources:
                            reference.sources.append(source)

    session.commit()
    print('regulondb', session.query(Interaction).count())