with open('PAO1/PSICQUIC/MPIDB.txt') as csvfile:
    reader = csv.DictReader(csvfile, delimiter='\t')
    for row in reader:
        interactors = []

        if ((row['Taxid interactor A'].split('|')[0] != 'taxid:208964(pseae)') |
                (row['Taxid interactor B'].split('|')[0] != 'taxid:208964(pseae)')): continue

        if (session.query(Interactor).filter(
                Interactor.id == row['#ID(s) interactor A'].split(':')[1]).first() != None):
            interactors.append(session.query(Interactor).filter(
                Interactor.id == row['#ID(s) interactor A'].split(':')[1]).one())
        elif row['#ID(s) interactor A'].split(':')[0] == 'uniprotkb':
            if (session.query(Protein).filter(
                    Protein.uniprotkb == row['#ID(s) interactor A'].split(':')[1]).first() != None):
                interactors.append(session.query(Protein).filter(
                    Protein.uniprotkb == row['#ID(s) interactor A'].split(':')[1]).one())
        if (session.query(Interactor).filter(
                Interactor.id == row['ID(s) interactor B'].split(':')[1]).first() != None):
            interactors.append(session.query(Interactor).filter(
                Interactor.id == row['ID(s) interactor B'].split(':')[1]).one())
        elif row['ID(s) interactor B'].split(':')[0] == 'uniprotkb':
            if (session.query(Protein).filter(
                    Protein.uniprotkb == row['ID(s) interactor B'].split(':')[1]).first() != None):
                interactors.append(session.query(Protein).filter(
                    Protein.uniprotkb == row['ID(s) interactor B'].split(':')[1]).one())

        if (len(interactors) != 2): continue

        homogenous = (interactors[0] == interactors[1])
        interaction = session.query(Interaction).filter((Interaction.interactors.contains(interactors[0])),
                                                        (Interaction.interactors.contains(interactors[1])),
                                                        (Interaction.homogenous == homogenous)).first()
        if (interaction == None):
            interaction = Interaction(strain='PAO1', type=(interactors[0].type + '-' + interactors[1].type),
                                      homogenous=homogenous, interactors=interactors)
            session.add(interaction), session.commit()

        reference = InteractionReference(interaction_type=row['Interaction type(s)'].split('(')[1][:-1],
                                         interaction_id=interaction.id,
                                         confidence_score=row['Confidence value(s)'],
                                         pmid=row['Publication Identifier(s)'].split('med:')[1][:8],
                                         detection_method=
                                         row['Interaction detection method(s)'].split('(')[1][:-1],
                                         author_last_name=row['Publication 1st author(s)'].split(' ')[0],
                                         publication_date=row['Publication 1st author(s)'].split('(')[1][:-1],
                                         publication_ref=row['Publication 1st author(s)'],
                                         source_db=row['Source database(s)'].split('(')[1][:-1],
                                         experimental_role_a=
                                         row['Experimental role(s) interactor A'].split('(')[1][:-1],
                                         experimental_role_b=
                                         row['Experimental role(s) interactor B'].split('(')[1][:-1],
                                         interactor_a_id=row['#ID(s) interactor A'].split(':')[1],
                                         interactor_b_id=row['ID(s) interactor B'].split(':')[1])
        session.add(reference)
        interaction.is_experimental = 1
    session.commit()
    print(session.query(Interaction).count())