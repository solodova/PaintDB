import csv
from DB_schema import Interactor, Protein, Interaction, InteractionReference, InteractionSource, InteractionXref
from DB_build import is_experimental_psimi

def parse_irefindex(session):
    parse_irefindex('PAO1/PSICQUIC/iRefIndex.txt', 'PAO1', 'taxid:208964(Pseudomonas aeruginosa PAO1)', session)
    parse_irefindex('PA14/PSICQUIC/iRefIndex.txt', 'PA14', 'taxid:208963(Pseudomonas aeruginosa UCBPP-PA14)', session)

def parse_irefindex(file, strain, taxid, session):
    with open(file) as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
            interactors = []

            if ((row['Taxid interactor A'].split('|')[0] != taxid) |
                    (row['Taxid interactor B'].split('|')[0] != taxid)): continue

            A_id = row['#ID(s) interactor A'].split(':')
            B_id = row['ID(s) interactor B'].split(':')

            if A_id[0] == 'uniprotkb':
                if session.query(Interactor).filter(Interactor.id == A_id[1]).first() is not None:
                    interactors.append(session.query(Interactor).filter(Interactor.id == A_id[1]).one())
                elif session.query(Protein).filter(Protein.uniprotkb == A_id[1]).first() is not None:
                    interactors.append(session.query(Protein).filter(Protein.uniprotkb == A_id[1]).one())
            elif A_id[0] == 'refseq':
                if session.query(Protein).filter(Protein.ncbi_acc == A_id[1]).first() is not None:
                    interactors.append(session.query(Protein).filter(Protein.ncbi_acc == A_id[1]).one())

            if B_id[0] == 'uniprotkb':
                if session.query(Interactor).filter(Interactor.id == B_id[1]).first() is not None:
                    interactors.append(session.query(Interactor).filter(Interactor.id == B_id[1]).one())
                elif session.query(Protein).filter(Protein.uniprotkb == B_id[1]).first() is not None:
                    interactors.append(session.query(Protein).filter(Protein.uniprotkb == B_id[1]).one())
            elif B_id[0] == 'refseq':
                if session.query(Protein).filter(Protein.ncbi_acc == B_id[1]).first() is not None:
                    interactors.append(session.query(Protein).filter(Protein.ncbi_acc == B_id[1]).one())

            if len(interactors) != 2: continue
            homogenous = (interactors[0] == interactors[1])

            interaction = session.query(Interaction).filter((Interaction.interactors.contains(interactors[0])),
                                                            (Interaction.interactors.contains(interactors[1])),
                                                            (Interaction.homogenous == homogenous)).first()
            if interaction is None:
                type = interactors[0].type + '-' + interactors[1].type
                interaction = Interaction(strain=strain, type=type, homogenous=homogenous, interactors=interactors)
                if row['Interaction detection method(s)'] != '-':
                    if is_experimental_psimi(row['Interaction detection method(s)'].split('MI:')[1][:4]):
                        interaction.is_experimental = 1
                    else:
                        interaction.is_experimental = 0
            else:
                if is_experimental_psimi(row['Interaction detection method(s)'].split('MI:')[1][:4]):
                    interaction.is_experimental = 1
                elif (row['Interaction detection method(s)'] == '-') and (interaction.is_experimental == 0):
                    interaction.is_experimental = None

            author, date, type= None, None, None
            pmids, detections = [None], [None]
            if row['Interaction detection method(s)'] != '-':
                del detections[0]
                for method in row['Interaction detection method(s)'].split('|'):
                    detections.append(method.split('(')[1][:-1])
            if (row['Interaction type(s)'] != '-'):
                type = row['Interaction type(s)'].split('(')[1][:-1]
            if (row['Publication 1st author(s)'] != '-'):
                author = row['Publication 1st author(s)'].split('-')[0][0].upper() + \
                         row['Publication 1st author(s)'].split('-')[0][1:]
                date = row['Publication 1st author(s)'].split('-')[1]
            if (row['Publication Identifier(s)'] != '-'):
                del pmids[0]
                for pmid in row['Publication Identifier(s)'].split('|'):
                    pmids.append(pmid.split('pubmed:')[1][:8])

            for pmid in pmids:
                for detection in detections:
                    reference = InteractionReference(interaction_id=interaction.id,
                                                     detection_method=detection,
                                                     author_ln=author,
                                                     pub_date=date,
                                                     pmid=pmid,
                                                     interaction_type=type,
                                                     source_db=row['Source database(s)'].split('(')[1][:-1],
                                                     confidence_score=row['Confidence value(s)'])
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
                                                             InteractionSource.data_source == 'iRefIndex').first()

            if source is None:
                source = InteractionSource(interaction_id=interaction.id, data_source='iRefIndex')
                session.add(source)
        print(session.query(Interaction).count())