import csv
from Schema1 import Interactor, OrthologPseudomonas, Interaction

inparalogs = {}

def parse_ortholuge(session):
    get_inparalogs('Orthologs/PAO1-PA14.csv')
    parse_orthologs('Orthologs/PAO1-PA14.csv', session)
    parse_ortholog_interactions(session)

def get_inparalogs(ortholog_file):
    with open(ortholog_file) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if (row['Strain 1 Inparalogs (Locus Tag/Name)'] != ''):
                strain1_inparalogs = row['Strain 1 Inparalogs (Locus Tag/Name)'].split(']')

                for inparalog in strain1_inparalogs:
                    if inparalog == '': continue
                    trimmed_inparalog = inparalog[:6]
                    if inparalog[0] == ';':
                        trimmed_inparalog = inparalog[1:7]

                    if trimmed_inparalog in inparalogs.keys():
                        inparalogs[trimmed_inparalog].append(row['Locus Tag (Strain 2)'])
                    else:
                        inparalogs[trimmed_inparalog] = [row['Locus Tag (Strain 2)']]

            if (row['Strain 2 Inparalogs (Locus Tag/Name)'] != ''):
                strain2_inparalogs = row['Strain 2 Inparalogs (Locus Tag/Name)'].split(']')

                for inparalog in strain2_inparalogs:
                    if inparalog == '': continue
                    trimmed_inparalog = inparalog[:9]
                    if inparalog[0] == ';':
                        trimmed_inparalog = inparalog[1:10]

                    if row['Locus Tag (Strain 1)'] in inparalogs.keys():
                        inparalogs[row['Locus Tag (Strain 1)']].append(trimmed_inparalog)
                    else:
                        inparalogs[row['Locus Tag (Strain 1)']] = [trimmed_inparalog]

def parse_orthologs(ortholog_file, session):
    with open(ortholog_file) as csvfile:
        reader2 = csv.DictReader(csvfile)
        for row in reader2:
            if (row['Locus Tag (Strain 1)']) in inparalogs.keys():
                if (row['Locus Tag (Strain 2)']) in inparalogs[row['Locus Tag (Strain 1)']]:
                    continue
            if (session.query(Interactor).filter(Interactor.id == row['Locus Tag (Strain 1)']).first() != None):
                if (session.query(Interactor).filter(Interactor.id == row['Locus Tag (Strain 2)']).first() != None):
                    ortholog1 = OrthologPseudomonas(protein_id=row['Locus Tag (Strain 1)'], strain_protein='PAO1',
                                                    ortholog_id=row['Locus Tag (Strain 2)'], strain_ortholog='PA14')
                    ortholog2 = OrthologPseudomonas(protein_id=row['Locus Tag (Strain 2)'], strain_protein='PA14',
                                                    ortholog_id=row['Locus Tag (Strain 1)'], strain_ortholog='PAO1')
                    session.add(ortholog1), session.add(ortholog2)
        session.commit()


def parse_ortholog_interactions(session):
    for interaction in session.query(Interaction).all():
        interactors, ortholog_interactors = [], [[], []]
        num = 0
        ortho = []
        for interactor in interaction.interactors:
            if interactor.type == 'protein':
                for ortholog in interactor.pseudomonas_orthologs:
                    if ortholog is not None:
                        ortholog_interactor = session.query(Interactor).filter(
                            Interactor.id == ortholog.ortholog_id).one()
                        ortholog_interactors[num].append(ortholog_interactor)
            else:
                ortholog_interactors[num].append(interactor)
            num += 1

        for interactor1 in ortholog_interactors[0]:
            for interactor2 in ortholog_interactors[1]:
                if (interactor1.type != 'metabolite') | (interactor2.type != 'metabolite'):
                    interactors.append([interactor1, interactor2])
        ortho.append(interactors)
        for interactor_pair in interactors:

            homogenous = (interactor_pair[0] == interactor_pair[1])
            ortholog_interaction = session.query(Interaction).filter(
                (Interaction.interactors.contains(interactor_pair[0])),
                (Interaction.interactors.contains(interactor_pair[1])),
                (Interaction.homogenous == homogenous)).first()

            if (ortholog_interaction != None):
                if (ortholog_interaction.ortholog_derived == None):
                    ortholog_interaction.ortholog_derived = 'confirmed from ' + interaction.strain
                elif ('from ' + interaction.strain) not in ortholog_interaction.ortholog_derived:
                    ortholog_interaction.ortholog_derived += ', confirmed from ' + interaction.strain
                session.commit()
            else:
                strain = None
                if interactor_pair[0].type == 'protein':
                    strain = interactor_pair[0].strain
                else:
                    strain = interactor_pair[1].strain
                ortholog_interaction = Interaction(strain=strain, type=interaction.type,
                                                   ortholog_derived='from ' + interaction.strain,
                                                   interactors=interactor_pair, homogenous=homogenous,
                                                   is_experimental=0)
                session.add(ortholog_interaction), session.commit()
    print(session.query(Interaction).count())