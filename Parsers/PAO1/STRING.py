with open('PAO1/STRING.txt') as csvfile:
    fieldnames = ['interactor_A', 'interactor_B', 'altID_A', 'altID_B', 'alias_A', 'alias_B', 'detection',
                  'publication', 'publication_ID', 'taxid_A', 'taxid_B', 'type', 'source_db', 'identifier',
                  'confidence']
    reader = csv.DictReader(csvfile, delimiter='\t', fieldnames=fieldnames)
    for row in reader:
        interactors = []
        locus_tag1 = row['interactor_A'].split('|')[0].split('.')[1]
        if (session.query(Interactor).filter(Interactor.id == locus_tag1).first() != None):
            interactors.append(session.query(Interactor).filter(Interactor.id == locus_tag1).one())
        locus_tag2 = row['interactor_B'].split('|')[0].split('.')[1]
        if (session.query(Interactor).filter(Interactor.id == locus_tag2).first() != None):
            interactors.append(session.query(Interactor).filter(Interactor.id == locus_tag2).one())

        if (len(interactors) != 2): continue
        homogenous = (interactors[0] == interactors[1])
        interaction = session.query(Interaction).filter((Interaction.interactors.contains(interactors[0])),
                                                        (Interaction.interactors.contains(interactors[1])),
                                                        (Interaction.homogenous == homogenous)).first()
        if (interaction == None):
            interaction = Interaction(strain='PAO1', type=(interactors[0].type + '-' + interactors[1].type),
                                      homogenous=homogenous, interactors=interactors)
            session.add(interaction), session.commit()

        detection = row['detection'].split('(')[1][:-1]
        confidence = row['confidence'].split(':')[1]
        type = find_type(row['type'].split(':')[2][:-1])
        source_db, author, date, ref, pmid = None, None, None, None, None

        if (row['source_db'] != '-'):
            source_db = row['source_db'].split('(')[1][:-1]
        if (row['publication'] != '-'):
            author = row['publication'].split(' ')[0]
            date = row['publication'].split('(')[1][:-1]
            ref = row['publication']
        if (row['publication_ID'] != '-'):
            pmid = row['publication_ID'].split(':')[1]
        reference = InteractionReference(interaction_type=type, interaction_id=interaction.id,
                                         detection_method=detection, confidence_score=confidence,
                                         source_db=source_db, pmid=pmid, author_last_name=author,
                                         publication_ref=ref, publication_date=date)
        session.add(reference)
        if (detection == 'experimental interaction detection'):
            interaction.is_experimental = 1
    session.commit()
    print(session.query(Interaction).count())