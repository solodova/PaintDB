import csv
from Schema1 import Interactor, Protein, Interaction, InteractionReference

def parse_regulatory_network(session):
    with open('regulatory_network.csv') as csvfile:
        reader = csv.DictReader(csvfile)
        string = []
        for row in reader:
            strains = row['Strain'].split(',')
            for strain in strains:
                if ((strain != 'PAO1') & (strain != 'PA14')): continue
                interactors = []
                if (session.query(Protein).filter((Protein.name == row['Regulator (TF or sigma)']),
                                                  (Protein.strain == strain)).first() is not None):
                    interactors.append(
                        session.query(Protein).filter((Protein.name == row['Regulator (TF or sigma)']),
                                                      (Protein.strain == strain)).one())
                elif (session.query(Interactor).filter(
                        Interactor.id == row['Regulator (TF or sigma)']).first() is not None):
                    interactors.append(session.query(Interactor).filter(
                        Interactor.id == row['Regulator (TF or sigma)']).one())
                if (session.query(Protein).filter((Protein.name == row['Target']),
                                                  (Protein.strain == strain)).first() is not None):
                    interactors.append(session.query(Protein).filter((Protein.name == row['Target']),
                                                                     (Protein.strain == strain)).one())
                elif (session.query(Interactor).filter(Interactor.id == row['Target']).first() is not None):
                    interactors.append(
                        session.query(Interactor).filter(Interactor.id == row['Target']).one())

                string.append(len(interactors))
                if (len(interactors) != 2): continue

                # string.append(row['Regulator (TF or sigma)'] + ':' + row['Target'] + '(' + strain + ')')
                homogenous = (interactors[0] == interactors[1])
                interaction = session.query(Interaction).filter((Interaction.interactors.contains(interactors[0])),
                                                                (Interaction.interactors.contains(interactors[1])),
                                                                (Interaction.homogenous == homogenous)).first()
                interactors[0].is_TF = 1

                if (interaction == None):
                    interaction = Interaction(strain=strain, type=interactors[0].type + ' - ' + interactors[1].type,
                                              homogenous=homogenous, interactors=interactors)
                    session.add(interaction), session.commit()

                reference = InteractionReference(interaction_id=interaction.id, pmid=row['pmid'],
                                                 interaction_type='TF/sigma-binding site (' + row['mode'] +
                                                                  'regulation)',
                                                 detection_method=row['evidence'], source_db=row['source_db'],
                                                 interaction_full_name=interactors[0].id + ' regulates(' +
                                                                       row['mode'] + ') ' + interactors[1].id)
                session.add(reference)
                interaction.is_experimental = 1
        session.commit()
        print(session.query(InteractionReference).count())