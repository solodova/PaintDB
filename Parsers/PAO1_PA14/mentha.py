import csv
from Schema1 import Interactor, Interaction, Protein, InteractionReference

def parse_mentha(session):
    parse_mentha('PAO1/PSICQUIC/mentha.txt', 'PAO1', 'taxid:208964(Pseudomonas aeruginosa PAO1)', session)
    parse_mentha('PA14/PSICQUIC/mentha.txt', 'PA14', 'taxid:208963(Pseudomonas aeruginosa UCBPP-PA14)', session)

def parse_mentha(file, strain, taxid, session):
    with open(file) as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
            interactors = []

            if ((row['Taxid interactor A'].split('|')[0] != taxid) |
                    (row['Taxid interactor B'].split('|')[0] != taxid)): continue

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
                interaction = Interaction(strain=strain, type=(interactors[0].type + '-' + interactors[1].type),
                                          homogenous=homogenous, interactors=interactors)
                interaction.interactors = interactors
                session.add(interaction), session.commit()

            reference = InteractionReference(interaction_type=row['Interaction type(s)'].split('(')[1][:-1],
                                             interaction_id=interaction.id,
                                             confidence_score=row['Confidence value(s)'],
                                             pmid=row['Publication Identifier(s)'].split('med:')[1][:8],
                                             detection_method=
                                             row['Interaction detection method(s)'].split('(')[1][:-1],
                                             source_db=row['Source database(s)'].split('(')[1][:-1])
            interaction.is_experimental = 1
            session.add(reference), session.commit()
        print(session.query(Interaction).count())