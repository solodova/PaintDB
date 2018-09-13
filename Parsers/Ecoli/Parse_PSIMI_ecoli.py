import csv
from Schema import Metabolite, Interaction, InteractionReference, OrthologEcoli, InteractionSource
from Parsing_Helpers import get_psimi_ref_list, is_experimental_interaction

# universal column names to be used for all file parsing
cols = ['interactor_A', 'interactor_B', 'altID_A', 'altID_B', 'alias_A', 'alias_B', 'detection', 'publication',
        'publication_ID', 'taxid_A', 'taxid_B', 'type', 'source_db', 'identifier', 'confidence']

# function to parse all Ecoli psimi formatted data files
def parse(session):
    parse_psimi(session, 'Data/Ecoli/PSICQUIC/BindingDB.txt', 'BindingDB(Ecoli)')
    parse_psimi(session, 'Data/Ecoli/PSICQUIC/EBI-GOA-nonIntAct.txt', 'EBI-GOA-nonIntAct(Ecoli)')
    parse_psimi(session, 'Data/Ecoli/PSICQUIC/IMEx.txt', 'IMEx(Ecoli)')
    parse_psimi(session, 'Data/Ecoli/PSICQUIC/IntAct.txt', 'IntAct(Ecoli)')
    parse_psimi(session, 'Data/Ecoli/PSICQUIC/iRefIndex.txt', 'iRefIndex(Ecoli)')
    parse_psimi(session, 'Data/Ecoli/PSICQUIC/mentha.txt', 'mentha(Ecoli)')
    parse_psimi(session, 'Data/Ecoli/PSICQUIC/MINT.txt', 'MINT(Ecoli)')
    parse_psimi(session, 'Data/Ecoli/PSICQUIC/MPIDB.txt', 'MPIDB(Ecoli)')
    parse_psimi(session, 'Data/Ecoli/PSICQUIC/UniProt.txt', 'UniProt(Ecoli)')
    parse_psimi(session, 'Data/Ecoli/DIP.txt', 'DIP(Ecoli)')
    parse_psimi(session, 'Data/Ecoli/IntAct.txt', 'IntAct(Ecoli)')

def parse_psimi(session, file, source):
    with open(file) as csvfile:
        reader = csv.DictReader(csvfile, fieldnames=cols, delimiter = '\t')

        # iterate through each interaction
        for row in reader:

            uniprot_A, refseq_A, orthologs_A , uniprot_B, refseq_B, orthologs_B= None, None, None, None, None, None
            # if one of the interactors is metabolite, save it's ids in pubchem and chebi
            pubchem, chebi = None, None
            # if one of the interactors is a metabolite, metabolite will be that metabolite and orthologs
            # will be set to the interaction's protein ortholog(s)
            metabolite_info, metabolite, orthologs = None, None, None

            # check if interactor A has uniprot or refseq id
            if 'uniprotkb' in row['interactor_A']:
                uniprot_A = row['interactor_A'].split('uniprotkb:')[1].split('|')[0]
            if 'refseq' in row['interactor_A']:
                refseq_A = row['interactor_A'].split('refseq:')[1].split('|')[0]

            # if uniprot id was found, look for orthologs matching that id
            if uniprot_A is not None:
                orthologs_A = session.query(OrthologEcoli).filter_by(ortholog_uniprot = uniprot_A).all()
            # if no orthologs were found but a refseq id was found, try to find ortholog based on refseq
            if (orthologs_A is None) and (refseq_A is not None):
                orthologs_A = session.query(OrthologEcoli).filter_by(ortholog_refseq = refseq_A).all()
            # if no orthologs were found for interactor A, but a uniprot or refseq does exist,
            # that means the ecoli interactor A is a protein without orthologs, so continue to next interaction
            if (orthologs_A is None) & ((uniprot_A is not None) | (refseq_A is not None)): continue

            # same as for interactor A above
            if 'uniprotkb' in row['interactor_B']:
                uniprot_B = row['interactor_B'].split('uniprotkb:')[1].split('|')[0]
            if 'refseq' in row['interactor_B']:
                refseq_B = row['interactor_B'].split('refseq:')[1].split('|')[0]

            if uniprot_B is not None:
                orthologs_B = session.query(OrthologEcoli).filter_by(ortholog_uniprot = uniprot_B).all()
            if (orthologs_B is None) and (refseq_B is not None):
                orthologs_B = session.query(OrthologEcoli).filter_by(ortholog_refseq = refseq_B).all()
            if (orthologs_B is None) & ((uniprot_B is not None) | (refseq_B is not None)): continue

            # if both orthologs_A and orthologs_B are None, then there are no protein interactors for this
            # interaction, so move on to the next interaction
            if (orthologs_A is None) and (orthologs_B is None): continue

            # if there were no orthologs for interactor A (and no refseq or uniprot was found),
            # search the file for pubchem or chebi ids for interactor A (as it may be a metabolite)
            if orthologs_A is None:
                if 'chebi' in row['interactor_A']:
                    chebi = row['interactor_A'].split('CHEBI:')[1].split('|')[0][:-1]
                if 'pubchem' in row['altID_A']:
                    pubchem = row['altID_A'].split('pubchem:')[1].split('|')[0]
                if (chebi is None) & ('chebi' in row['altID_A']):
                    chebi = row['altID_A'].split('CHEBI:')[1].split('|')[0][:-1]
                # if no metabolite ids were found in the interaction row, then move on to the next interaction
                # because no interactor_A was identified
                if (chebi is None) & (pubchem is None): continue
                # if a pubchem or chebi id was found, then this interaction will be a p-m interaction, so
                # set the protein interactors(orthologs) to orthologs_B
                orthologs = orthologs_B
            # other case where orthologs_B were not identified so need to check if interactor B has metabolite ids
            elif orthologs_B is None:
                if 'chebi' in row['interactor_B']:
                    chebi = row['interactor_B'].split('CHEBI:')[1].split('|')[0][:-1]
                if 'pubchem' in row['altID_B']:
                    pubchem = row['altID_B'].split('pubchem:')[1].split('|')[0]
                if (chebi is None) & ('chebi' in row['altID_B']):
                    chebi = row['altID_B'].split('CHEBI:')[1].split('|')[0][:-1]
                if (chebi is None) & (pubchem is None): continue
                orthologs = orthologs_A

            # if one of the interactors was identified to be a metabolite, search for the metabolite and set metabolite
            # variable to that value. if the metabolite doesnt exist create it
            # Note: if this point was reached, it means one of the interactors had protein orthologs,
            # so we can safely create a new metabolite knowing it will have a protein interaction partner
            if (chebi is not None) | (pubchem is not None):
                id = None
                # preferentially set id for new metabolites to be chebi
                if chebi is not None:
                    id = chebi
                    metabolite = session.query(Metabolite).filter_by(chebi = chebi).first()
                # if no metabolite with chebi was found, but pubchem id exists, try to find
                # metabolite with that pubchem
                if (metabolite is None) & (pubchem is not None):
                    id = pubchem
                    metabolite = session.query(Metabolite).filter_by(pubchem = pubchem).first()
                # if no metabolite was found with pubchem or chebi id, create new metabolite
                if metabolite is None:
                    metabolite = Metabolite(id = id, chebi=chebi, pubchem=pubchem)
                    session.add(metabolite)
                # if a metabolite was found, update its chebi and pubchem if it has none
                else:
                    if metabolite.pubchem is None:
                        metabolite.pubchem = pubchem
                    if metabolite.chebi is None:
                        metabolite.chebi = chebi

            # list of interactor pairs for interaction
            interactors = []
            # if no metabolite was found for interaction, it is a p-p interaction, so iterate through
            # orthologs to create interactor pairs
            if metabolite is None:
                for ortholog_A in orthologs_A:
                    for ortholog_B in orthologs_B:
                        if (ortholog_A is not None) and (ortholog_B is not None):
                            # only add the interactor pair if the protein strains match
                            if ortholog_A.strain_protein == ortholog_B.strain_protein:
                                interactors.append([[ortholog_A.protein, ortholog_A.ortholog_id],
                                                    [ortholog_B.protein, ortholog_B.ortholog_id]])
            else:
                # if a metabolite was found, add pairs of all orthologs with metabolite to interactor pairs
                for ortholog in orthologs:
                    interactors.append([[metabolite, metabolite.id], [ortholog.protein, ortholog.ortholog_id]])

            # for each interactor pair, create interaction if it doesnt exist, otherwise update attributes
            for interactor_pair in interactors:
                homogenous = (interactor_pair[0][0] == interactor_pair[1][0])
                interaction = session.query(Interaction).filter(Interaction.interactors.contains(interactor_pair[0][0]),
                                                                Interaction.interactors.contains(interactor_pair[1][0]),
                                                                Interaction.homogenous == homogenous).first()
                if interaction is None:
                    # since one of the interactors may be a metabolite, set strain to match strain of protein
                    strain = None
                    if interactor_pair[0][0].type == 'p':
                        strain = interactor_pair[0][0].strain
                    else:
                        strain = interactor_pair[1][0].strain
                    # if interaction did not exist, set it to Ecoli ortholog derived
                    interaction = Interaction(strain=strain,
                                              interactors=[interactor_pair[0][0], interactor_pair[1][0]],
                                              type=(interactor_pair[0][0].type + '-' + interactor_pair[1][0].type),
                                              ortholog_derived='Ecoli')
                    session.add(interaction), session.commit()

                ref_parameter_list = get_psimi_ref_list(row)

                # in case the interaction already existed, make sure interactor_a and interactor_b variables for
                # new interaction reference match up with the first and second interactors of the existing
                # interaction (so it's easy to see which Pseudomonas interactor matches up with which Ecoli
                # ortholog)
                interactor_a, interactor_b = None, None
                if interaction.interactors[0] == interactor_pair[0][0]:
                    interactor_a = interactor_pair[0][1]
                    interactor_b = interactor_pair[1][1]
                else:
                    interactor_b = interactor_pair[0][1]
                    interactor_a = interactor_pair[1][1]

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
                        detection_method=ref[0], author_ln = ref[1], pub_date=ref[2], pmid=ref[3],
                        interaction_type=ref[4], source_db=ref[5], confidence=ref[6],
                        interactor_a=interactor_a, interactor_b=interactor_b).first()
                    # if nref doesn't exist, create and add it to the interaction's reference list,
                    # and add the source to the reference's sources
                    if nref is None:
                        nref = InteractionReference(
                            detection_method=ref[0], author_ln = ref[1], pub_date=ref[2], pmid=ref[3],
                            interaction_type=ref[4], source_db=ref[5], confidence=ref[6],
                            interactor_a=interactor_a, interactor_b=interactor_b)
                        interaction.references.append(nref)
                        nref.sources.append(nsource)
                    # if nref does exist, add the interaction and source to it's attributes if they aren't added
                    else:
                        if interaction not in nref.interactions:
                            nref.interactions.append(interaction)
                        if nsource not in nref.sources:
                            nref.sources.append(nsource)

    session.commit()
    print(source, session.query(Interaction).count())