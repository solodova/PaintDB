from Schema import Metabolite, Protein, Interaction, InteractionSource
import csv
from sqlalchemy import func, or_


def get_db_stats_interactors(strain, session):
    # new_file = open('PaIntDB_PAO1_stats', mode='x')
    header = ['Statistic', 'Value']
    file_name = 'PaIntDB_' + strain + '_Interactors.csv'
    file_writer = csv.DictWriter(open(file_name, mode='x', newline=''), fieldnames=header)
    file_writer.writeheader()
    file_writer.writerow(
        {'Statistic': 'Num proteins(p)',
         'Value': session.query(Protein).filter(Protein.strain == strain,
                                                or_(Protein.uniprotkb == None, Protein.uniprotkb != 'pc')).count()})
    file_writer.writerow(
        {'Statistic': 'Num proteins(pc)',
         'Value': session.query(Protein).filter(Protein.strain == strain, Protein.uniprotkb == 'pc').count()})
    file_writer.writerow(
        {'Statistic': 'Num p with uniprotkb',
         'Value': session.query(Protein).filter(Protein.strain == strain, Protein.uniprotkb != 'pc',
                                                Protein.uniprotkb != None).count()})
    file_writer.writerow(
        {'Statistic': 'Num p with refseq',
         'Value': session.query(Protein).filter(Protein.strain == strain,
                                                or_(Protein.uniprotkb != 'pc', Protein.uniprotkb == None),
                                                Protein.ncbi_acc != None).count()})
    file_writer.writerow(
        {'Statistic': 'Num p with name',
         'Value': session.query(Protein).filter(Protein.strain == strain,
                                                or_(Protein.uniprotkb != 'pc', Protein.uniprotkb == None),
                                                Protein.name != None).count()})
    file_writer.writerow(
        {'Statistic': 'Num p with product name',
         'Value': session.query(Protein).filter(Protein.strain == strain,
                                                or_(Protein.uniprotkb != 'pc', Protein.uniprotkb == None),
                                                Protein.product_name != None).count()})
    # file_writer.writerow({'Statistic': 'Num p with >0 localization',
    #                       'Value': session.query(Protein).filter(Protein.strain == 'PAO1').join(Protein.localizations).
    #                       group_by(Protein).having(func.count(Protein.localizations) > 0).count()}))
    # file_writer.writerow({'Statistic': 'Num p with ==1 localization',
    #                       'Value': session.query(Protein).filter(Protein.strain == 'PAO1').join(Protein.localizations).
    #                      group_by(Protein).having(func.count(Protein.localizations) == 1).count()}))
    # file_writer.writerow({'Statistic': 'Num p with >1 localization',
    #                       'Value': session.query(Protein).filter(Protein.strain == 'PAO1').join(Protein.localizations).
    #                      group_by(Protein).having(func.count(Protein.localizations) > 1).count()}))
    file_writer.writerow({'Statistic': 'Num p with >0 ontology',
                          'Value': session.query(Protein).filter(Protein.strain == strain).join(Protein.ontologies).
                         group_by(Protein).having(func.count(Protein.ontologies) > 0).count()})
    file_writer.writerow({'Statistic': 'Num p with ==1 ontology',
                          'Value': session.query(Protein).filter(Protein.strain == strain).join(Protein.ontologies).
                         group_by(Protein).having(func.count(Protein.ontologies) == 1).count()})
    file_writer.writerow({'Statistic': 'Num p with >1 ontology',
                          'Value': session.query(Protein).filter(Protein.strain == strain).join(Protein.ontologies).
                         group_by(Protein).having(func.count(Protein.ontologies) > 1).count()})
    file_writer.writerow({'Statistic': 'Proteins with >0 ortholog in Pseudomonas (p)',
                          'Value': session.query(Protein).filter(Protein.strain == strain).join(
                              Protein.pseudomonas_orthologs).
                         group_by(Protein).having(func.count(Protein.pseudomonas_orthologs) > 0).count()})
    file_writer.writerow({'Statistic': 'Proteins with ==1 ortholog in Pseudomonas (p)',
                          'Value': session.query(Protein).filter(Protein.strain == strain).join(
                              Protein.pseudomonas_orthologs).
                         group_by(Protein).having(func.count(Protein.pseudomonas_orthologs) == 1).count()})
    file_writer.writerow({'Statistic': 'Proteins with >1 ortholog in Pseudomonas (p)',
                          'Value': session.query(Protein).filter(Protein.strain == strain).join(
                              Protein.pseudomonas_orthologs).
                         group_by(Protein).having(func.count(Protein.pseudomonas_orthologs) > 1).count()})
    file_writer.writerow({'Statistic': 'Proteins with >0 ortholog in Ecoli (p)',
                          'Value': session.query(Protein).filter(Protein.strain == strain).join(
                              Protein.ecoli_orthologs).
                         group_by(Protein).having(func.count(Protein.ecoli_orthologs) > 0).count()})
    file_writer.writerow({'Statistic': 'Proteins with ==1 ortholog in Ecoli (p)',
                          'Value': session.query(Protein).filter(Protein.strain == strain).join(
                              Protein.ecoli_orthologs).
                         group_by(Protein).having(func.count(Protein.ecoli_orthologs) == 1).count()})
    file_writer.writerow({'Statistic': 'Proteins with >1 ortholog in Ecoli (p)',
                          'Value': session.query(Protein).filter(Protein.strain == strain).join(
                              Protein.ecoli_orthologs).
                         group_by(Protein).having(func.count(Protein.ecoli_orthologs) > 1).count()})
    # file_writer.writerow({'Statistic': 'Num p, pc with interactions',
    #                       'Value': session.query(Protein).filter(Protein.strain == 'PAO1').join(Protein.interactions).
    #                      group_by(Protein).having(func.count(Protein.interactions) > 0).count()}))
    # file_writer.writerow({'Statistic': 'Num p, pc with p-p interactions',
    #                       'Value': session.query(Protein).filter(Protein.strain == 'PAO1').join(Protein.interactions).
    #                       filter(Interaction.type == 'p-p').
    #                      group_by(Protein).having(func.count(Protein.interactions)> 0).count()}))
    # file_writer.writerow({'Statistic': 'Num p, pc with p-m interactions',
    #                       'Value': session.query(Protein).filter(Protein.strain == 'PAO1').join(Protein.interactions).
    #                       filter(Interaction.type.in_(['p-m', 'm-p'])).
    #                      group_by(Protein).having(func.count(Protein.interactions) > 0).count()}))

def get_db_stats_interactions(strain, session):
    header = ['Statistic', 'Total', 'Total P-P', 'Total P-M', 'Experimental', 'Experimental P-P',
              'Experimental P-M', 'Non-experimental', 'Non-experimental P-P', 'Non-experimental P-M',
              'Unknown detection', 'Unknown detection P-P', 'Unknown detection P-M']
    file_name = 'PaIntDB_' + strain + '_Interaction.csv'
    file_writer = csv.DictWriter(open(file_name, mode='x', newline=''), fieldnames=header)
    file_writer.writeheader()
    file_writer.writerow(
        {'Statistic': 'Num interactions',
         'Total':
             session.query(Interaction).filter_by(strain=strain).count(),
         'Total P-P':
             session.query(Interaction).filter_by(strain=strain, type='p-p').count(),
         'Total P-M':
             session.query(Interaction).filter(Interaction.strain == strain, Interaction.type != 'p-p').count(),
         'Experimental':
             session.query(Interaction).filter_by(strain=strain).join('sources').
                 filter(InteractionSource.is_experimental == 1).count(),
         'Experimental P-P':
             session.query(Interaction).filter_by(strain=strain, type='p-p').join('sources').
                 filter(InteractionSource.is_experimental == 1).count(),
         'Experimental P-M':
             session.query(Interaction).filter(Interaction.strain == strain, Interaction.type != 'p-p').
                 join('sources').filter(InteractionSource.is_experimental == 1).count(),
         'Non-experimental':
             session.query(Interaction).filter_by(strain=strain).join('sources').
                 filter(InteractionSource.is_experimental == 0).count(),
         'Non-experimental P-P':
             session.query(Interaction).filter_by(strain=strain, type='p-p').join('sources').
                 filter(InteractionSource.is_experimental == 0).count(),
         'Non-experimental P-M':
             session.query(Interaction).filter(Interaction.strain == strain, Interaction.type != 'p-p').
                 join('sources').filter(InteractionSource.is_experimental == 0).count(),
         'Unknown detection':
             session.query(Interaction).filter_by(strain=strain).join('sources').
                 filter(InteractionSource.is_experimental == 2).count(),
         'Unknown detection P-P':
             session.query(Interaction).filter_by(strain=strain, type='p-p').join('sources').
                 filter(InteractionSource.is_experimental == 2).count(),
         'Unknown detection P-M':
             session.query(Interaction).filter(Interaction.strain == strain, Interaction.type != 'p-p').
                 join('sources').filter(InteractionSource.is_experimental == 2).count()})

    print('done')
    sources_Ecoli = ['EcoCyc', 'RegulonDB(Ecoli)', 'IMEx(Ecoli)', 'BindingDB(Ecoli)', 'EBI-GOA-nonIntAct(Ecoli)',
                     'IntAct(Ecoli)', 'iRefIndex(Ecoli)', 'mentha(Ecoli)', 'MINT(Ecoli)', 'MPIDB(Ecoli)',
                     'UniProt(Ecoli)', 'DIP(Ecoli)', 'KEGG(Ecoli)']

    sources_PA14 = ['IMEx(PA14)', 'IntAct(PA14)', 'iRefIndex(PA14)', 'mentha(PA14)', 'MINT(PA14)', 'KEGG(PA14)',
                    'Galan-Vasquez(PA14)']
    sources_PAO1 = ['Geoff', 'XLinkDB', 'Zhang', 'ADIPInteractomes(PAO1)', 'IMEx(PAO1)', 'IntAct(PAO1)',
                    'iRefIndex(PAO1)', 'mentha(PAO1)', 'MINT(PAO1)', 'Galan-Vasquez(PAO1)', 'KEGG(PAO1)']

    for source in [['PAO1', sources_PAO1], ['PA14', sources_PA14], ['Ecoli', sources_Ecoli]]:
        file_writer.writerow(
            {'Statistic': 'Num interactions from ' + source[0] + ' sources',
             'Total':
                 session.query(Interaction).filter_by(strain=strain).join('sources').group_by(Interaction).
                     filter(InteractionSource.data_source.in_(source[1])).count(),
             'Total P-P':
                 session.query(Interaction).filter_by(strain=strain, type='p-p').join('sources').
                     group_by(Interaction).filter(InteractionSource.data_source.in_(source[1])).count(),
             'Total P-M':
                 session.query(Interaction).filter(Interaction.strain == strain, Interaction.type != 'p-p').
                     join('sources').group_by(Interaction).filter(
                     InteractionSource.data_source.in_(source[1])).count(),
             'Experimental':
                 session.query(Interaction).filter_by(strain=strain).join('sources').group_by(Interaction).
                     filter(InteractionSource.is_experimental == 1,
                            InteractionSource.data_source.in_(source[1])).count(),
             'Experimental P-P':
                 session.query(Interaction).filter_by(strain=strain, type='p-p').join('sources').
                     group_by(Interaction).filter(InteractionSource.is_experimental == 1,
                                                  InteractionSource.data_source.in_(source[1])).count(),
             'Experimental P-M':
                 session.query(Interaction).filter(Interaction.strain == strain, Interaction.type != 'p-p').
                     join('sources').group_by(Interaction).
                     filter(InteractionSource.is_experimental == 1,
                            InteractionSource.data_source.in_(source[1])).count(),
             'Non-experimental':
                 session.query(Interaction).filter_by(strain=strain).join('sources').group_by(Interaction).
                     filter(InteractionSource.is_experimental == 0,
                            InteractionSource.data_source.in_(source[1])).count(),
             'Non-experimental P-P':
                 session.query(Interaction).filter_by(strain=strain, type='p-p').join('sources').
                     group_by(Interaction).filter(InteractionSource.is_experimental == 0,
                                                  InteractionSource.data_source.in_(source[1])).count(),
             'Non-experimental P-M':
                 session.query(Interaction).filter(Interaction.strain == strain, Interaction.type != 'p-p').
                     join('sources').group_by(Interaction).
                     filter(InteractionSource.is_experimental == 0,
                            InteractionSource.data_source.in_(source[1])).count(),
             'Unknown detection':
                 session.query(Interaction).filter_by(strain=strain).join('sources').group_by(Interaction).
                     filter(InteractionSource.is_experimental == 2,
                            InteractionSource.data_source.in_(source[1])).count(),
             'Unknown detection P-P':
                 session.query(Interaction).filter_by(strain=strain, type='p-p').join('sources').
                     group_by(Interaction).filter(InteractionSource.is_experimental == 2,
                                                  InteractionSource.data_source.in_(source[1])).count(),
             'Unknown detection P-M':
                 session.query(Interaction).filter(Interaction.strain == strain, Interaction.type != 'p-p').
                     join('sources').group_by(Interaction).
                     filter(InteractionSource.is_experimental == 2,
                            InteractionSource.data_source.in_(source[1])).count()})
    print('done1')
    all_sources = sources_PAO1 + sources_PA14 + sources_Ecoli
    for source in all_sources:
        file_writer.writerow(
            {'Statistic': 'Num interactions from ' + source,
             'Total':
                 session.query(InteractionSource).filter_by(data_source=source).join('interactions').
                     filter(Interaction.strain == strain).count(),
             'Total P-P':
                 session.query(InteractionSource).filter_by(data_source=source).join('interactions').
                     filter(Interaction.strain == strain, Interaction.type == 'p-p').count(),
             'Total P-M':
                 session.query(InteractionSource).filter_by(data_source=source).join('interactions').
                     filter(Interaction.strain == strain, Interaction.type != 'p-p').count(),
             'Experimental':
                 session.query(InteractionSource).filter_by(data_source=source, is_experimental=1).
                     join('interactions').filter(Interaction.strain == strain).count(),
             'Experimental P-P':
                 session.query(InteractionSource).filter_by(data_source=source, is_experimental=1).
                     join('interactions').filter(Interaction.strain == strain,
                                                 Interaction.type == 'p-p').count(),
             'Experimental P-M':
                 session.query(InteractionSource).filter_by(data_source=source, is_experimental=1).
                     join('interactions').filter(Interaction.strain == strain,
                                                 Interaction.type != 'p-p').count(),
             'Non-experimental':
                 session.query(InteractionSource).filter_by(data_source=source, is_experimental=0).
                     join('interactions').filter(Interaction.strain == strain).count(),
             'Non-experimental P-P':
                 session.query(InteractionSource).filter_by(data_source=source, is_experimental=0).
                     join('interactions').filter(Interaction.strain == strain,
                                                 Interaction.type == 'p-p').count(),
             'Non-experimental P-M':
                 session.query(InteractionSource).filter_by(data_source=source, is_experimental=0).
                     join('interactions').filter(Interaction.strain == strain,
                                                 Interaction.type != 'p-p').count(),
             'Unknown detection':
                 session.query(InteractionSource).filter_by(data_source=source, is_experimental=2).
                     join('interactions').filter(Interaction.strain == strain).count(),
             'Unknown detection P-P':
                 session.query(InteractionSource).filter_by(data_source=source, is_experimental=2).
                     join('interactions').filter(Interaction.strain == strain,
                                                 Interaction.type == 'p-p').count(),
             'Unknown detection P-M':
                 session.query(InteractionSource).filter_by(data_source=source, is_experimental=2).
                     join('interactions').filter(Interaction.strain == strain,
                                                 Interaction.type != 'p-p').count()})
    print('done2')
    p_sources = sources_PAO1
    p_ortholog_sources = sources_PA14
    if strain == 'PA14':
        p_sources = sources_PA14
    p_ortholog_sources = sources_PAO1

    # for source in [['Pseudomonas', p_ortholog_sources, p_sources],
    #                ['Ecoli', sources_Ecoli, sources_PAO1 + sources_PA14]]:
    #     file_writer.writerow(
    #         {'Statistic': 'Num interactions in ' + strain + ' from ' + source[0] + ' orthology',
    #          'Total':
    #              session.query(Interaction).filter_by(strain=strain).join('sources').group_by(Interaction).
    #                  filter(Interaction.sources.any(InteractionSource.data_source.in_(source[1]))).
    #                  filter(~(Interaction.sources.any(InteractionSource.data_source.in_(source[2])))).count(),
    #          'Total P-P':
    #              session.query(Interaction).filter_by(strain=strain, type='p-p').
    #         join('sources').group_by(Interaction).
    #         filter(Interaction.sources.any(InteractionSource.data_source.in_(source[1]))).
    #         filter(~(Interaction.sources.any(InteractionSource.data_source.in_(source[2])))).count(),
    #          'Total P-M':
    #              session.query(Interaction).filter(Interaction.strain == strain, Interaction.type != 'p-p').
    #         join('sources').group_by(Interaction).
    #         filter(Interaction.sources.any(InteractionSource.data_source.in_(source[1]))).
    #         filter(~(Interaction.sources.any(InteractionSource.data_source.in_(source[2])))).count(),
    #          'Experimental':
    #              session.query(Interaction).filter_by(strain=strain).join('sources').group_by(Interaction).
    #         filter(InteractionSource.is_experimental == 1).
    #         filter(Interaction.sources.any(InteractionSource.data_source.in_(source[1]))).
    #         filter(~(Interaction.sources.any(InteractionSource.data_source.in_(source[2])))).count(),
    #          'Experimental P-P':
    #              session.query(Interaction).filter_by(strain=strain, type='p-p').
    #         join('sources').group_by(Interaction).filter(InteractionSource.is_experimental == 1).
    #         filter(Interaction.sources.any(InteractionSource.data_source.in_(source[1]))).
    #         filter(~(Interaction.sources.any(InteractionSource.data_source.in_(source[2])))).count(),
    #          'Experimental P-M':
    #              session.query(Interaction).filter(Interaction.strain == strain, Interaction.type != 'p-p').
    #         join('sources').group_by(Interaction).filter(InteractionSource.is_experimental == 1).
    #         filter(Interaction.sources.any(InteractionSource.data_source.in_(source[1]))).
    #         filter(~(Interaction.sources.any(InteractionSource.data_source.in_(source[2])))).count(),
    #          'Non-experimental':
    #              session.query(Interaction).filter_by(strain=strain).join('sources').group_by(Interaction).
    #         filter(InteractionSource.is_experimental == 0).
    #         filter(Interaction.sources.any(InteractionSource.data_source.in_(source[1]))).
    #         filter(~(Interaction.sources.any(InteractionSource.data_source.in_(source[2])))).count(),
    #          'Non-experimental P-P':
    #              session.query(Interaction).filter_by(strain=strain, type='p-p').
    #         join('sources').group_by(Interaction).filter(InteractionSource.is_experimental == 0).
    #         filter(Interaction.sources.any(InteractionSource.data_source.in_(source[1]))).
    #         filter(~(Interaction.sources.any(InteractionSource.data_source.in_(source[2])))).count(),
    #          'Non-experimental P-M':
    #              session.query(Interaction).filter(Interaction.strain == strain, Interaction.type != 'p-p').
    #         join('sources').group_by(Interaction).filter(InteractionSource.is_experimental == 0).
    #         filter(Interaction.sources.any(InteractionSource.data_source.in_(source[1]))).
    #         filter(~(Interaction.sources.any(InteractionSource.data_source.in_(source[2])))).count(),
    #          'Unknown detection':
    #              session.query(Interaction).filter_by(strain=strain).join('sources').group_by(Interaction).
    #         filter(InteractionSource.is_experimental == 2).
    #         filter(Interaction.sources.any(InteractionSource.data_source.in_(source[1]))).
    #         filter(~(Interaction.sources.any(InteractionSource.data_source.in_(source[2])))).count(),
    #          'Unknown detection P-P':
    #              session.query(Interaction).filter_by(strain=strain, type='p-p').
    #         join('sources').group_by(Interaction).filter(InteractionSource.is_experimental == 2).
    #         filter(Interaction.sources.any(InteractionSource.data_source.in_(source[1]))).
    #         filter(~(Interaction.sources.any(InteractionSource.data_source.in_(source[2])))).count(),
    #          'Unknown detection P-M':
    #              session.query(Interaction).filter(Interaction.strain == strain, Interaction.type != 'p-p').
    #         join('sources').group_by(Interaction).filter(InteractionSource.is_experimental == 2).
    #         filter(Interaction.sources.any(InteractionSource.data_source.in_(source[1]))).
    #         filter(~(Interaction.sources.any(InteractionSource.data_source.in_(source[2])))).count()})


    # for source in [['Pseudomonas', p_ortholog_sources, p_sources],
    #                ['Ecoli', sources_Ecoli, sources_PAO1 + sources_PA14]]:
    #     file_writer.writerow(
    #         {'Statistic': 'Num interactions in ' + strain + ' from ' + source[0] + ' orthology',
    #          'Total':
    #              session.query(Interaction)
    #                  .filter_by(strain=strain)
    #                  .join('sources')
    #                  .group_by(Interaction)
    #                  .filter(Interaction.sources.any(InteractionSource.data_source.in_(source[1]))).
    #                  filter(~(Interaction.sources.any(InteractionSource.data_source.in_(source[2]))))
    #                  .count(),
    #          'Total P-P':'-',
    #          'Total P-M':'-',
    #          'Experimental':'-',
    #          'Experimental P-P':'-',
    #          'Experimental P-M':'-',
    #          'Non-experimental':'-',
    #          'Non-experimental P-P':'-',
    #          'Non-experimental P-M':'-',
    #          'Unknown detection':'-',
    #          'Unknown detection P-P':'-',
    #          'Unknown detection P-M': '-'})
