import csv
from DB_schema import Interactor, Interaction, Protein, InteractionReference, InteractionXref, InteractionSource
from DB_build import is_experimental_psimi

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

            A_id = row['#ID(s) interactor A'].split(':')[1]
            B_id = row['ID(s) interactor B'].split(':')[1]

            if session.query(Interactor).filter(Interactor.id == A_id).first() is not None:
                interactors.append(session.query(Interactor).filter(Interactor.id == A_id).one())
            elif session.query(Protein).filter(Protein.uniprotkb == A_id).first() is not None:
                interactors.append(session.query(Protein).filter(Protein.uniprotkb == A_id).one())

            if session.query(Interactor).filter(Interactor.id == B_id).first() is not None:
                interactors.append(session.query(Interactor).filter(Interactor.id == B_id).one())
            elif session.query(Protein).filter(Protein.uniprotkb == B_id).first() is not None:
                interactors.append(session.query(Protein).filter(Protein.uniprotkb == B_id).one())

            if len(interactors) != 2: continue
            homogenous = (interactors[0] == interactors[1])

            interaction = session.query(Interaction).filter(Interaction.interactors.contains(interactors[0]),
                                                            Interaction.interactors.contains(interactors[1]),
                                                            Interaction.homogenous == homogenous).first()
            if interaction is None:
                type=(interactors[0].type + '-' + interactors[1].type)
                interaction = Interaction(strain=strain, type=type, homogenous=homogenous, interactors=interactors)
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
                                             pmid=row['Publication Identifier(s)'].split('pubmed:')[1][:8],
                                             interaction_type=row['Interaction type(s)'].split('(')[1][:-1],
                                             source_db=row['Source database(s)'].split('(')[1][:-1],
                                             confidence_score=row['Confidence value(s)'])
            session.add(reference)

            xref_field = row['Interaction identifier(s)'].split(':')
            xref = session.query(InteractionXref).filter(InteractionXref.accession == xref_field[1],
                                                         InteractionXref.interaction_id == interaction.id).first()

            if xref is None:
                xref = InteractionXref(interaction_id=interaction.id, accession=xref_field[1],
                                       data_source=xref_field[0])
                session.add(xref)

            source = session.query(InteractionSource).filter(InteractionSource.interaction_id == interaction.id,
                                                             InteractionSource.data_source == 'mentha').first()

            if source is None:
                source = InteractionSource(interaction_id=interaction.id, data_source='mentha')
                session.add(source)
        print(session.query(Interaction).count())