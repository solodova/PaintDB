import csv

from Schema1 import Interactor, Protein, Interaction, InteractionReference
from Parsers.Parser import is_experimental_psimi

def parse_ADIPInteractomes(session):
    with open('PAO1/PSICQUIC/ADIPInteractomes.txt') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
            interactors = []

            if ((row['Taxid interactor A'] != 'taxid:208964(Pseudomonas aeruginosa PAO1)') |
                    (row['Taxid interactor B'] != 'taxid:208964(Pseudomonas aeruginosa PAO1)')): continue

            if (session.query(Interactor).filter(
                    Interactor.id == row['#ID(s) interactor A'].split(':')[1]).first() != None):
                interactors.append(session.query(Interactor).filter(
                    Interactor.id == row['#ID(s) interactor A'].split(':')[1]).one())
            elif (session.query(Protein).filter(
                    Protein.uniprotkb == row['#ID(s) interactor A'].split(':')[1]).first() != None):
                interactors.append(session.query(Protein).filter(
                    Protein.uniprotkb == row['#ID(s) interactor A'].split(':')[1]).one())
            if (session.query(Interactor).filter(
                    Interactor.id == row['ID(s) interactor B'].split(':')[1]).first() != None):
                interactors.append(session.query(Interactor).filter(
                    Interactor.id == row['ID(s) interactor B'].split(':')[1]).one())
            elif (session.query(Protein).filter(
                    Protein.uniprotkb == row['ID(s) interactor B'].split(':')[1]).first() != None):
                interactors.append(session.query(Protein).filter(
                    Protein.uniprotkb == row['ID(s) interactor B'].split(':')[1]).one())

            if (len(interactors) != 2): continue
            homogenous = (interactors[0] == interactors[1])
            interaction = session.query(Interaction).filter((Interaction.interactors.contains(interactors[0])),
                                                            (Interaction.interactors.contains(interactors[1])),
                                                            (Interaction.homogenous == homogenous)).first()
            if (interaction == None):
                interaction = Interaction(strain='PAO1', homogenous=homogenous, interactors=interactors)
                session.add(interaction), session.commit()

            reference = InteractionReference(interaction_type=row['Interaction type(s)'],
                                             interaction_id=interaction.id,
                                             confidence_score=row['Confidence value(s)'],
                                             pmid=row['Publication Identifier(s)'].split('med:')[1][:8],
                                             detection_method=
                                             row['Interaction detection method(s)'].split('(')[1][:-1],
                                             author_last_name=row['Publication 1st author(s)'].split(',')[0],
                                             publication_date=
                                             row['Publication 1st author(s)'].split('(')[1][:-1],
                                             publication_ref=row['Publication 1st author(s)'])
            session.add(reference)
            #if is_experimental_psimi()interaction.is_experimental = 1
            # xref = session.query(InteractionXref).filter(
            #     InteractionXref.accession == row['Interaction identifier(s)'].split(':')[1]).first()
            # if (xref == None):
            #     xref = InteractionXref(accession=row['Interaction identifier(s)'].split(':')[1],
            #                            data_source=row['Source database(s)'].split('(')[1][:-1],
            #                            interaction_id=interaction.id)
            #     session.add(xref)
        session.commit()
        print(session.query(Interaction).count())
