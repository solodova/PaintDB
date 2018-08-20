import csv
from Schema1 import Interactor, Interaction, OrthologEcoli, InteractionReference

def parse_Ecoli_DIP(session):

    with open('Ecoli/DIP.txt') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
            interactors = []

            ids_A = row['ID interactor A'].split('|')
            ids_B = row['ID interactor B'].split('|')
            refseq_A, uniprotkb_A, refseq_B, uniprotkb_B = '', '', '', ''
            for id in ids_A:
                fields = id.split(':')
                if (fields[0] == 'refseq'):
                    refseq_A = fields[1]
                elif (fields[0] == 'uniprotkb'):
                    uniprotkb_A = fields[1]
            for id in ids_B:
                fields = id.split(':')
                if (fields[0] == 'refseq'):
                    refseq_B = fields[1]
                elif (fields[0] == 'uniprotkb'):
                    uniprotkb_B = fields[1]

            orthologs_A, orthologs_B = [], []
            if (uniprotkb_A != ''):
                orthologs_A = session.query(OrthologEcoli).filter(OrthologEcoli.ortholog_uniprot == uniprotkb_A).all()
            if ((len(orthologs_A) == 0) & (refseq_A != '')):
                orthologs_A = session.query(OrthologEcoli).filter(OrthologEcoli.ortholog_refseq == refseq_A).all()
            if (uniprotkb_B != ''):
                orthologs_B = session.query(OrthologEcoli).filter(OrthologEcoli.ortholog_uniprot == uniprotkb_B).all()
            if ((len(orthologs_B) == 0) & (refseq_B != '')):
                orthologs_B = session.query(OrthologEcoli).filter(OrthologEcoli.ortholog_refseq == refseq_B).all()

            for ortholog_A in orthologs_A:
                if ortholog_A != None:
                    interactor_A = session.query(Interactor).filter(Interactor.id == ortholog_A.protein_id).one()
                    for ortholog_B in orthologs_B:
                        if ortholog_B != None:
                            interactor_B = session.query(Interactor).filter(
                                Interactor.id == ortholog_B.protein_id).one()
                            if (interactor_A.strain == interactor_B.strain):
                                interactors.append([interactor_A, interactor_B])

            for interactor_pair in interactors:
                homogenous = (interactor_pair[0] == interactor_pair[1])
                interaction = session.query(Interaction).filter(
                    (Interaction.interactors.contains(interactor_pair[0])),
                    (Interaction.interactors.contains(interactor_pair[1])),
                    (Interaction.homogenous == homogenous)).first()
                if (interaction != None):
                    if (interaction.ortholog_derived == None):
                        interaction.ortholog_derived = 'confirmed from E.coli'
                    elif ('from E. coli' not in interaction.ortholog_derived):
                        interaction.ortholog_derived += ', confirmed from E. coli'
                    session.commit()
                if (interaction == None):
                    interaction = Interaction(strain=interactor_pair[0].strain, interactors=interactor_pair,
                                              type=(interactor_pair[0].type + '-' + interactor_pair[1].type),
                                              is_experimental=0, ortholog_derived='from E. coli')
                    session.add(interaction), session.commit()
                detections, pmids, types, list = [], [], [], []
                if (row['Interaction detection method(s)'] != '-'):
                    detections = row['Interaction detection method(s)'].split('|')
                    list.append(detections)
                if (row['Publication Identifier(s)'] != '-'):
                    pmids = row['Publication Identifier(s)'].split('|')
                    list.append(pmids)
                if (row['Interaction type(s)']):
                    types = row['Interaction type(s)'].split('|')
                    list.append(types)

                if (len(list) != 0) & (all((len(item) == len(list[0])) for item in list)):
                    for num in range(0, len(list[0])):
                        type, pmid, detection = None, None, None
                        for item in list:
                            if (item == types):
                                type = types[num].split('(')[1][:-1]
                            if (item == pmids):
                                pmid = pmids[num * 2].split('med:')[1][:8]
                            if (item == detections):
                                detection = detections[num].split('(')[1][:-1]
                        # there are more than one pmid sometimes
                        reference = InteractionReference(interaction_type=type, interaction_id=interaction.id,
                                                         pmid=pmid, detection_method=detection,
                                                         source_db=row['Source database(s)'].split('(')[1][:-1])
                        session.add(reference)
        session.commit()
        print(session.query(Interaction).count())