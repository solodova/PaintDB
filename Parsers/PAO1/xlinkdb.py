with open('PAO1/xlinkdb.txt') as csvfile:
    reader = csv.DictReader(csvfile, delimiter='\t')
    for row in reader:
        interactors = []
        if (session.query(Interactor).filter(Interactor.id == row['proA']).first() != None):
            interactors.append(session.query(Interactor).filter(Interactor.id == row['proA']).one())
        elif (session.query(Protein).filter(Protein.uniprotkb == row['proA']).first() != None):
            interactors.append(session.query(Protein).filter(Protein.uniprotkb == row['proA']).one())
        if (session.query(Interactor).filter(Interactor.id == row['proB']).first() != None):
            interactors.append(session.query(Interactor).filter(Interactor.id == row['proB']).one())
        elif (session.query(Protein).filter(Protein.uniprotkb == row['proB']).first() != None):
            interactors.append(session.query(Protein).filter(Protein.uniprotkb == row['proB']).one())

        if (len(interactors) != 2): continue
        homogenous = (interactors[0] == interactors[1])
        interaction = session.query(Interaction).filter((Interaction.interactors.contains(interactors[0])),
                                                        (Interaction.interactors.contains(interactors[1])),
                                                        (Interaction.homogenous == homogenous)).first()
        if (interaction == None):
            interaction = Interaction(strain='PAO1', type=(interactors[0].type + '-' + interactors[1].type),
                                      homogenous=homogenous, interactors=interactors)
            session.add(interaction), session.commit()
        reference = InteractionReference(interaction_type='physical association', interaction_id=interaction.id,
                                         pmid='25800553', source_db='XLinkDB',
                                         detection_method='chemical cross-linking mass spectrometry')
        session.add(reference)
        interaction.is_experimental = 1
    session.commit()
    print(session.query(Interaction).count())