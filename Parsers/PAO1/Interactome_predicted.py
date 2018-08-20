with open('PAO1/Interactome_predicted.csv') as csvfile:
    reader = csv.DictReader(csvfile)

    for row in reader:
        if (float(row['Confidence']) < 0.9): continue
        interactors = []
        if (session.query(Interactor).filter(Interactor.id == row['Protein1']).first() != None):
            interactors.append(session.query(Interactor).filter(Interactor.id == row['Protein1']).one())

        if (session.query(Interactor).filter(Interactor.id == row['Protein2']).first() != None):
            interactors.append(session.query(Interactor).filter(Interactor.id == row['Protein2']).one())

        if (len(interactors) != 2): continue
        homogenous = (interactors[0] == interactors[1])
        interaction = session.query(Interaction).filter((Interaction.interactors.contains(interactors[0])),
                                                        (Interaction.interactors.contains(interactors[1])),
                                                        (Interaction.homogenous == homogenous)).first()
        if (interaction == None):
            interaction = Interaction(strain='PAO1', type=(interactors[0].type + '-' + interactors[1].type),
                                      homogenous=homogenous, is_experimental=0, interactors=interactors)
            session.add(interaction), session.commit()

            reference = InteractionReference(interaction_type='predicted', interaction_id=interaction.id,
                                             confidence_score=row['Confidence'], pmid='22848443',
                                             detection_method='computational prediction')
            session.add(reference)
    session.commit()
    print(session.query(Interaction).count())