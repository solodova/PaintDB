import csv
from Schema1 import Interactor, Protein, Interaction, InteractionReference, ProteinXref

def parse_iRefIndex(session):
    parse_iRefIndex('PAO1/PSICQUIC/iRefIndex.txt', 'PAO1', 'taxid:208964(Pseudomonas aeruginosa PAO1)', session)
    parse_iRefIndex('PA14/PSICQUIC/iRefIndex.txt', 'PA14', 'taxid:208963(Pseudomonas aeruginosa UCBPP-PA14)', session)

def parse_iRefIndex(file, strain, taxid, session):
    with open(file) as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
            interactors = []

            if ((row['Taxid interactor A'].split('|')[0] != taxid) |
                    (row['Taxid interactor B'].split('|')[0] != taxid)): continue

            if (row['#ID(s) interactor A'].split(':')[0] == 'uniprotkb'):
                if (session.query(Interactor).filter(
                        Interactor.id == row['#ID(s) interactor A'].split(':')[1]).first() != None):
                    interactors.append(session.query(Interactor).filter(
                        Interactor.id == row['#ID(s) interactor A'].split(':')[1]).one())
                else:
                    if (session.query(Protein).filter(
                            Protein.uniprotkb == row['#ID(s) interactor A'].split(':')[1]).first() != None):
                        interactors.append(session.query(Protein).filter(
                            Protein.uniprotkb == row['#ID(s) interactor A'].split(':')[1]).one())
            elif (row['#ID(s) interactor A'].split(':')[0] == 'refseq'):
                if (session.query(ProteinXref).filter(ProteinXref.accession ==
                                                      row['#ID(s) interactor A'].split(':')[1]).first() != None):
                    interactors.append(session.query(ProteinXref).filter(ProteinXref.accession ==
                                                                         row['#ID(s) interactor A'].split(
                                                                             ':')[1]).one().protein)
            if (row['ID(s) interactor B'].split(':')[0] == 'uniprotkb'):
                if (session.query(Interactor).filter(
                        Interactor.id == row['ID(s) interactor B'].split(':')[1]).first() != None):
                    interactors.append(session.query(Interactor).filter(
                        Interactor.id == row['ID(s) interactor B'].split(':')[1]).one())
                elif (session.query(Protein).filter(
                        Protein.uniprotkb == row['ID(s) interactor B'].split(':')[1]).first() != None):
                    interactors.append(session.query(Protein).filter(
                        Protein.uniprotkb == row['ID(s) interactor B'].split(':')[1]).one())
            elif (row['ID(s) interactor B'].split(':')[0] == 'refseq'):
                if (session.query(ProteinXref).filter(ProteinXref.accession ==
                                                      row['ID(s) interactor B'].split(':')[1]).first() != None):
                    interactors.append(session.query(ProteinXref).filter(ProteinXref.accession ==
                                                                         row['ID(s) interactor B'].split(
                                                                             ':')[1]).one().protein)

            if (len(interactors) != 2): continue

            homogenous = (interactors[0] == interactors[1])
            interaction = session.query(Interaction).filter((Interaction.interactors.contains(interactors[0])),
                                                            (Interaction.interactors.contains(interactors[1])),
                                                            (Interaction.homogenous == homogenous)).first()
            if (interaction == None):
                interaction = Interaction(strain=strain, type=(interactors[0].type + '-' + interactors[1].type),
                                          homogenous=homogenous, interactors=interactors)
                session.add(interaction), session.commit()

            author, date, ref, type, pmid, detection = None, None, None, None, None, None
            if (row['Publication 1st author(s)'] != '-'):
                author = row['Publication 1st author(s)'].split('-')[0]
                date = row['Publication 1st author(s)'].split('-')[1]
                ref = row['Publication 1st author(s)']
            if (row['Interaction type(s)'] != '-'):
                type = row['Interaction type(s)'].split('(')[1][:-1]
            if (row['Publication Identifier(s)'] != '-'):
                pmid = row['Publication Identifier(s)'].split('med:')[1][:8]
            if (row['Interaction detection method(s)'] != '-'):
                detection = row['Interaction detection method(s)'].split('(')[1][:-1]
            # there are more than one pmid sometimes
            reference = InteractionReference(interaction_type=type, interaction_id=interaction.id, pmid=pmid,
                                             confidence_score=row['Confidence value(s)'], author_last_name=author,
                                             detection_method=detection, publication_date=date, publication_ref=ref,
                                             source_db=row['Source database(s)'].split('(')[1][:-1],
                                             experimental_role_a=
                                             row['Experimental role(s) interactor A'].split('(')[1][:-1],
                                             experimental_role_b=
                                             row['Experimental role(s) interactor B'].split('(')[1][:-1],
                                             interactor_a_id=row['#ID(s) interactor A'].split(':')[1],
                                             interactor_b_id=row['ID(s) interactor B'].split(':')[1])
            interaction.is_experimental = 1
            session.add(reference), session.commit()
        print(session.query(Interaction).count())