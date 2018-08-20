with open('PAO1/GeoffWinsor.csv') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        interactors = []
        if (session.query(Interactor).filter(Interactor.id == row['locus_tag']).first() != None):
            interactors.append(session.query(Interactor).filter(Interactor.id == row['locus_tag']).one())
        row = next(reader)
        if (session.query(Interactor).filter(Interactor.id == row['locus_tag']).first() != None):
            interactors.append(session.query(Interactor).filter(Interactor.id == row['locus_tag']).one())

        if (len(interactors) != 2): continue
        homogenous = (interactors[0] == interactors[1])
        interaction = session.query(Interaction).filter((Interaction.interactors.contains(interactors[0])),
                                                        (Interaction.interactors.contains(interactors[1])),
                                                        (Interaction.homogenous == homogenous)).first()
        if (interaction == None):
            interaction = Interaction(strain='PAO1', type=(interactors[0].type + '-' + interactors[1].type),
                                      homogenous=homogenous, interactors=interactors)
            session.add(interaction), session.commit()

        reference = InteractionReference(interaction_type=row['type'], interaction_id=interaction.id,
                                         pmid=row['pmid'], interaction_full_name=row['full_name'],
                                         detection_method=row['experimental_type'])
        interaction.is_experimental = 1
        session.add(reference)
    session.commit()
    print(session.query(Interaction).count())