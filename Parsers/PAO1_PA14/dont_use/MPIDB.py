import csv
from Schema import Interactor, Protein, Interaction, InteractionSource, InteractionReference, InteractionXref
from Main import is_experimental_psimi

def parse_mpidb(session):
    with open('PAO1/PSICQUIC/MPIDB.txt') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
            interactors = []

            if (row['Taxid interactor A'].split('|')[0] != 'taxid:208964(pseae)') |\
                    (row['Taxid interactor B'].split('|')[0] != 'taxid:208964(pseae)'): continue

            A_id = row['#ID(s) interactor A'].split(':')[1]
            B_id = row['ID(s) interactor B'].split(':')[1]

            if session.query(Interactor).filter(Interactor.id == A_id).first() is not None:
                interactors.append(session.query(Interactor).filter(Interactor.id ==A_id).one())
            elif session.query(Protein).filter(Protein.uniprotkb == A_id).first() is not None:
                interactors.append(session.query(Protein).filter(Protein.uniprotkb == A_id).one())

            if session.query(Interactor).filter(Interactor.id == B_id).first() is not None:
                interactors.append(session.query(Interactor).filter(Interactor.id == B_id).one())
            elif session.query(Protein).filter(Protein.uniprotkb == B_id).first() is not None:
                interactors.append(session.query(Protein).filter(Protein.uniprotkb == B_id).one())

            if len(interactors) != 2: continue
            homogenous = (interactors[0] == interactors[1])

            interaction = session.query(Interaction).filter((Interaction.interactors.contains(interactors[0])),
                                                            (Interaction.interactors.contains(interactors[1])),
                                                            (Interaction.homogenous == homogenous)).first()
            if interaction is None:
                type = interactors[0].type + '-' + interactors[1].type
                interaction = Interaction(strain='PAO1', type=type, homogenous=homogenous, interactors=interactors)
                if is_experimental_psimi(row['Interaction detection method(s)'].split('MI:')[1][:4]):
                    interaction.is_experimental = 1
                else:
                    interaction.is_experimental = 0
                session.add(interaction), session.commit()
            else:
                if is_experimental_psimi(row['Interaction detection method(s)'].split('MI:')[1][:4]):
                    interaction.is_experimental = 1

            reference = InteractionReference(interaction_id=interaction.id,
                                             detection_method=row['Interaction detection method(s)'].split('(')[1][:-1],
                                             author_ln=row['Publication 1st author(s)'].split(' ')[0],
                                             pub_date=row['Publication 1st author(s)'].split('(')[1][:-1],
                                             pmid=row['Publication Identifier(s)'].split('pubmed:')[1][:8],
                                             confidence=row['Confidence value(s)'],
                                             interaction_type=row['Interaction type(s)'].split('(')[1][:-1],
                                             source_db=row['Source database(s)'])
            session.add(reference)

            for xref in row['Interaction identifier(s)'].split('|'):
                xref_field = xref.split(':')
                xref = session.query(InteractionXref).filter(InteractionXref.accession == xref_field[1],
                                                             InteractionXref.interaction_id == interaction.id).first()

                if xref is None:
                    xref = InteractionXref(interaction_id=interaction.id, accession=xref_field[1],
                                           data_source=xref_field[0])
                    session.add(xref)

            source = session.query(InteractionSource).filter(InteractionSource.interaction_id == interaction.id,
                                                             InteractionSource.data_source == 'MPIDB').first()

            if source is None:
                source = InteractionSource(interaction_id=interaction.id, data_source='MPIDB')
                session.add(source)
        session.commit()
        print(session.query(Interaction).count())