import csv
from DB_schema import Interactor, Protein, Interaction, InteractionReference, InteractionSource, InteractionXref
from DB_parsing_helpers import get_psimi_ref_list, is_experimental_interaction

cols = ['interactor_A', 'interactor_B', 'altID_A', 'altID_B', 'alias_A', 'alias_B', 'detection', 'publication',
        'publication_ID', 'taxid_A', 'taxid_B', 'type', 'source_db', 'identifier', 'confidence']

def parse(session):
    parse_psimi('Data/PAO1/PSICQUIC/ADIPInteractomes.txt', 'PAO1', 'ADIPInteractomes(PAO1)', session)
    parse_psimi('Data/PAO1/PSICQUIC/IMEx.txt', 'PAO1', 'IMEx(PAO1)', session)
    parse_psimi('Data/PAO1/PSICQUIC/IntAct.txt', 'PAO1', 'IntAct(PAO1)', session)
    parse_psimi('Data/PAO1/PSICQUIC/iRefIndex.txt', 'PAO1', 'iRefIndex(PAO1)', session)
    parse_psimi('Data/PAO1/PSICQUIC/mentha.txt', 'PAO1', 'mentha(PAO1)', session)
    parse_psimi('Data/PAO1/PSICQUIC/MINT.txt', 'PAO1', 'MINT(PAO1)', session)
    parse_psimi('Data/PAO1/PSICQUIC/MPIDB.txt', 'PAO1', 'MPIDB(PAO1)', session)
    parse_psimi('Data/PAO1/IntAct.txt', 'PAO1', 'IntAct(PAO1)', session)

    parse_psimi('Data/PA14/PSICQUIC/IMEx.txt', 'PA14', 'IMEx(PA14)', session)
    parse_psimi('Data/PA14/PSICQUIC/iRefIndex.txt', 'PA14', 'iRefIndex(PA14)', session)
    parse_psimi('Data/PA14/PSICQUIC/mentha.txt', 'PA14', 'mentha(PA14)', session)
    parse_psimi('Data/PA14/IntAct.txt', 'PA14', 'IntAct(PA14)', session)

def parse_psimi(file, strain, source, session):
    with open(file) as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t', fieldnames=cols)
        next(reader)
        for row in reader:
            uniprot_A, refseq_A, interactor_A, uniprot_B, refseq_B, interactor_B = None, None, None, None, None, None

            # check if interactor A has uniprot or refseq id, store these values
            if 'uniprotkb' in row['interactor_A']:
                uniprot_A = row['interactor_A'].split('uniprotkb:')[1].split('|')[0]
            if 'refseq' in row['interactor_A']:
                refseq_A = row['interactor_A'].split('refseq:')[1].split('|')[0]

            # if a uniprot id was found, try to find the interactor in the database
            if uniprot_A is not None:
                # check if there is a protein-complex with this uniprot id
                interactor_A = session.query(Interactor).get(uniprot_A)
                # if no protein complex, check for protein matching the uniprot id
                if interactor_A is None:
                    interactor_A = session.query(Protein).filter_by(uniprotkb = uniprot_A).first()
            # if no interactor A was found but there was also a refseq id, try to find the protein based on
            # it's refseq
            if (interactor_A is None) and (refseq_A is not None):
                interactor_A = session.query(Protein).filter_by(ncbi_acc = refseq_A).first()
            # if no interactor A was found, move on to next interaction
            if interactor_A is None: continue

            # same as for A above
            if 'uniprotkb' in row['interactor_B']:
                uniprot_B = row['interactor_B'].split('uniprotkb:')[1].split('|')[0]
            if 'refseq' in row['interactor_B']:
                refseq_B = row['interactor_B'].split('refseq:')[1].split('|')[0]

            if uniprot_B is not None:
                interactor_B = session.query(Interactor).get(uniprot_B)
                if interactor_B is None:
                    interactor_B = session.query(Protein).filter_by(uniprotkb = uniprot_B).first()
            if (interactor_B is None) and (refseq_B is not None):
                interactor_B = session.query(Protein).filter_by(ncbi_acc = refseq_B).first()

            if interactor_B is None: continue


            homogenous = (interactor_A == interactor_B)
            interaction = session.query(Interaction).filter(Interaction.interactors.contains(interactor_A),
                                                            Interaction.interactors.contains(interactor_B),
                                                            Interaction.homogenous == homogenous).first()
            # if no interaction was found with the interactors, create a new interaction
            if interaction is None:
                interaction = Interaction(strain=strain, type='p-p', homogenous=homogenous,
                                          interactors=[interactor_A, interactor_B])
                session.add(interaction), session.commit()

            ref_parameter_list = get_psimi_ref_list(row)

            is_experimental = is_experimental_interaction(row)

            # check to see if source exists
            nsource = session.query(InteractionSource).filter_by(
                data_source=source, is_experimental=is_experimental).first()
            # if source doesn't exist, create and add it to the interaction's sources
            if nsource is None:
                nsource = InteractionSource(data_source=source, is_experimental=is_experimental)
                interaction.sources.append(nsource)
            # if the source does exist, add it to the interaction's sources if it isn't already
            elif nsource not in interaction.sources:
                interaction.sources.append(nsource)

            # go through each reference in the ref_parameter list, search for it, and if it doesnt exist create it
            for ref in ref_parameter_list:
                nref = session.query(InteractionReference).filter_by(
                    detection_method=ref[0], author_ln=ref[1], pub_date=ref[2], pmid=ref[3],
                    interaction_type=ref[4], source_db=ref[5], confidence=ref[6],
                    interactor_a=None, interactor_b=None).first()
                # if nref doesn't exist, create and add it to the interaction's reference list,
                # and add the source to the reference's sources
                if nref is None:
                    nref = InteractionReference(
                        detection_method=ref[0], author_ln=ref[1], pub_date=ref[2], pmid=ref[3],
                        interaction_type=ref[4], source_db=ref[5], confidence=ref[6])
                    interaction.references.append(nref)
                    nref.sources.append(nsource)
                # if nref does exist, add the interaction and source to it's attributes if they aren't added
                else:
                    if interaction not in nref.interactions:
                        nref.interactions.append(interaction)
                    if nsource not in nref.sources:
                        nref.sources.append(nsource)

            #collect all the cross references for the interaction
            for xref in row['identifier'].split('|'):
                xref_field = xref.split(':')
                # check if the cross reference exists for this interaction, if it doesnt create it
                xref = session.query(InteractionXref).filter_by(accession = xref_field[1],
                                                                interaction_id = interaction.id).first()

                if xref is None:
                    xref = InteractionXref(interaction_id=interaction.id, accession=xref_field[1],
                                           data_source=xref_field[0])
                    session.add(xref)

        session.commit()
    print(source, session.query(Interaction).count())