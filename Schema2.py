def parse_kegg(org_id, strain):
    # get pathways for organism specified by org_id
    pathways = kegg_list(database='pathway', org=org_id).read().split('path:')
    path_names, path_ids = [], []

    for path in pathways:
        if path != '':
            path_names.append(path.split('\t')[1].split(' -')[0])
            path_ids.append(path[:8])

    mas_compound_ids = []
    for path in path_ids:
        # get kgml representation of path
        kgml_path = read(kegg_get(path, option='kgml'))
        path_name = kgml_path._getname()
        # dictionary of compounds in current path (node_id: kegg_id)
        #   compound._getid() returns node id (only relevant in context of current path)
        #   compound._getname() returns kegg id (relevant in overall KEGG DB
        compound_ids = {}
        for compound in kgml_path.compounds:
            compound_ids[compound._getid()] = compound._getname()[-6:]
        mas_compound_ids.append(compound_ids)

        # go through each relation in path
        for relation in kgml_path.relations:
            relation_type = relation.element.attrib['type']

            if (relation_type == 'maplink'): continue
            # relation._getentry1/2() returns  protein id (locus) or compound id (KEGG id)
            entries = [relation._getentry1()._getname(), relation._getentry2()._getname()]
            if (entries[0] == 'undefined') | (entries[1] == 'undefined'): continue
            interactors = [[], []]
            new_metabolites = [[], []]
            # go through each entry in the relation, find/create interactors
            for num in range(0, 2):
                # each entry may contain >1 id; go through all of them
                for id in entries[num].split(' '):
                    if (id == ''): continue
                    # if interactor is not protein or compound, continue
                    if (id.split(':')[0] != org_id) & (id.split(':')[1] not in compounds): continue

                    # check if interactor (protein or metabolite) already exists
                    if (session.query(Interactor).filter(Interactor.id == id.split(':')[1]).first() is not None):
                        interactors[num].append(
                            session.query(Interactor).filter(Interactor.id == id.split(':')[1]).one())
                    # if it doesnt exist, it's not a valid protein, so check if it is a valid compound
                    elif (id.split(':')[1] in compounds):
                        # if it is a valid compound, create new metabolite
                        new_metabolites[num].append(id.split(':')[1])
                    # if parsing E. coli path, add all orthologs to interactor list
                    elif (org_id == 'eco'):
                        for ortholog in (session.query(OrthologEcoli).filter(
                                (OrthologEcoli.ortholog_id == id.split(':')[1]),
                                (OrthologEcoli.strain_protein == strain)).all()):
                            interactors[num].append(ortholog.protein)

            # create list of interactor pairs from two separate lists (interactors[0], interactors[1])
            interactor_pairs = []
            for interactor1 in interactors[0]:
                for interactor2 in interactors[1]:
                    if (interactor1.type != 'metabolite') | (interactor2.type != 'metabolite'):
                        interactor_pairs.append([interactor1, interactor2])

            for interactor1 in interactors[0]:
                for interactor2 in new_metabolites[1]:
                    if (interactor1.type != 'metabolite'):
                        new_metabolite = Metabolite(id=interactor2, type='metabolite',
                                                name=compounds[interactor2][0],
                                                pubchem_id=compounds[interactor2][1], kegg_id=interactor2)
                        session.add(new_metabolite), session.commit()
                        interactor_pairs.append([interactor1, new_metabolite])

            for interactor1 in interactors[1]:
                for interactor2 in new_metabolites[0]:
                    if (interactor1.type != 'metabolite'):
                        new_metabolite = Metabolite(id=interactor2, type='metabolite',
                                                name=compounds[interactor2][0],
                                                pubchem_id=compounds[interactor2][1], kegg_id=interactor2)
                        session.add(new_metabolite), session.commit()
                        interactor_pairs.append([interactor1, new_metabolite])

            if (len(interactor_pairs) == 0): continue

            # get all intermediates in reaction of type compound
            intermeds = []
            for subtype in relation.element.iter(tag='subtype'):
                if 'compound' in subtype.attrib:
                    compound_node_id = subtype.attrib['compound']
                    if (compound_node_id == None): continue
                    if (int(compound_node_id) not in compound_ids): continue
                    # if compound id is valid, either add existing matching metabolite or create new one and add
                    compound_id = compound_ids[int(compound_node_id)]
                    metabolite = session.query(Metabolite).filter(Metabolite.kegg_id == compound_id).first()
                    if metabolite == None:
                        metabolite = Metabolite(kegg_id=compound_id, pubchem_id=compounds[compound_id][1],
                                                name=compounds[compound_id][0], type='metabolite', id=compound_id)
                        session.add(metabolite)
                        session.commit()
                    else:
                        metabolite = session.query(Metabolite).filter(Metabolite.kegg_id == compound_id).one()
                    intermeds.append(metabolite)

            # add protein - intermediate interactor pairs
            for interactor_list in interactors:
                for interactor in interactor_list:
                    if (interactor.type != 'metabolite'):
                        for intermed in intermeds:
                            interactor_pairs.append([interactor, intermed])

            for interactor_pair in interactor_pairs:
                homogenous = (interactor_pair[0] == interactor_pair[1])
                interaction = session.query(Interaction).filter(
                    (Interaction.interactors.contains(interactor_pair[0])),
                    (Interaction.interactors.contains(interactor_pair[1])),
                    (Interaction.homogenous == homogenous)).first()

                if (interaction == None):
                    interaction = Interaction(type=interactor_pair[0].type + '-' + interactor_pair[1].type,
                                              strain=strain, is_experimental=0, homogenous=homogenous,
                                              interactors=interactor_pair)
                    if org_id == 'eco':
                        interaction.ortholog_derived = 'from E. coli'
                    session.add(interaction), session.commit()
                else:
                    interaction = session.query(Interaction).filter(
                        (Interaction.interactors.contains(interactor_pair[0])),
                        (Interaction.interactors.contains(interactor_pair[1])),
                        (Interaction.homogenous == homogenous)).one()
                    if org_id == 'eco':
                        if (interaction.ortholog_derived == None):
                            interaction.ortholog_derived = 'from E. coli'
                        else:
                            if ('from E.coli' not in interaction.ortholog_derived):
                                interaction.ortholog_derived += ', confirmed from E.coli'

                if (session.query(InteractionXref).filter((InteractionXref.interaction_id == interaction.id),
                                                          (InteractionXref.data_source == 'KEGG')).first() == None):
                    xref = InteractionXref(interaction_id=interaction.id, data_source='KEGG')
                    session.add(xref), session.commit()

    print(session.query(Interaction).count())
