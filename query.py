import os
from flask import Flask, Blueprint, flash, g, redirect, render_template, request, url_for, send_file
from .helpers import removeFiles
from werkzeug.utils import secure_filename
from flask import current_app as app
from .DB_schema import Interactor, Interaction, InteractionSource
import re, csv
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

engine = create_engine('sqlite:////Users/olga/Desktop/PaIntDB.db', echo=True)
Session = sessionmaker(bind=engine)

ALLOWED_EXTENSIONS = set(['txt', 'csv', 'tsv'])
bp = Blueprint('query', __name__, url_prefix='/query')

psimi_fields = ['ID(s) interactor A', 'ID(s) interactor B', 'Alt. ID(s) interactor A', 'Alt. ID(s) interactor B',
                'Alias(es) interactor A', 'Alias(es) interactor B',	'Interaction detection method(s)',
                'Publication 1st author(s)', 'Publication identifier(s)', 'Taxid interactor A', 'Taxid interactor B',
                'Interaction type(s)', 'Source database(s)', 'Interaction identifier(s)', 'Confidence value(s)',
                'Annotation(s) interactor A', 'Annotation(s) interactor B']

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

@bp.route('upload', methods=['GET', 'POST'])
def upload():
    if request.method == 'POST':
        # if no files were added before submit, flash error
        if 'proteinFile' not in request.files:
            flash('Error: No file selected.')
            return redirect(request.url)

        protein_file = request.files['proteinFile']
        if protein_file.filename == '':
            flash('Error: No file selected.')
            return redirect(request.url)
        if not allowed_file(protein_file.filename):
            flash('Error: Incompatible file type for ' + protein_file.filename)
            return redirect(request.url)

        if protein_file is not None:
            pfile = secure_filename(('proteins.' + protein_file.filename.rsplit('.', 1)[1].lower()))
            #pfile = 'proteins.' + protein_file.filename.rsplit('.', 1)[1].lower()
            parts = re.split(r'/', pfile)
            protein_file.save(os.path.join(app.config['UPLOAD_FOLDER'], *parts))
            flash(protein_file.filename + ' was successfully uploaded! Click below to move on to the next step.')
        return redirect(url_for('query.uploaded'))
    #removeFiles(app)
    return render_template('query/upload.html')

@bp.route('uploaded', methods=['GET', 'POST'])
def uploaded():
    if request.method =='POST':
        #removeFiles(app)
        return redirect(url_for('query.upload'))
    return render_template('query/uploaded.html')

@bp.route('filter', methods=['GET', 'POST'])
def filter():
    if request.method == 'POST':
        filters = {'strain': [], 'type': [], 'Ecoli_sources': [], 'PAO1_sources': [], 'PA14_sources': [],
                   'verification': []}
        sources_Ecoli = ['EcoCyc', 'RegulonDB(Ecoli)', 'IMEx(Ecoli)', 'BindingDB(Ecoli)', 'EBI-GOA-nonIntAct(Ecoli)',
                         'IntAct(Ecoli)', 'iRefIndex(Ecoli)', 'mentha(Ecoli)', 'MINT(Ecoli)', 'MPIDB(Ecoli)',
                         'UniProt(Ecoli)', 'DIP(Ecoli)', 'KEGG(Ecoli)']
        sources_PAO1 = ['Geoff', 'XLinkDB', 'Zhang', 'ADIPInteractomes(PAO1)', 'IMEx(PAO1)', 'IntAct(PAO1)',
                        'iRefIndex(PAO1)', 'mentha(PAO1)', 'MINT(PAO1)', 'Galan-Vasquez(PAO1)', 'KEGG(PAO1)']
        sources_PA14 = ['IMEx(PA14)', 'IntAct(PA14)', 'iRefIndex(PA14)', 'mentha(PA14)', 'KEGG(PA14)',
                        'Galan-Vasquez(PA14)']
        filters_all ={'strain': ['PAO1', 'PA14'], 'type': ['p-p', 'p-m', 'p-bs'],
                      'Ecoli_sources': sources_Ecoli, 'PAO1_sources': sources_PAO1, 'PA14_sources': sources_PA14,
                      'verification': ['0', '1', '2']}
        sources = []
        tfbs_sources = ['RegulonDB(Ecoli)', 'Galan-Vasquez(PAO1)']

        for filter in filters:
            if filter in request.form:
                filters[filter] = request.form.getlist(filter)
                if 'None' in filters[filter]:
                    filters[filter] = ['None']
                elif 'All' in filters[filter]:
                    filters[filter] = filters_all[filter]

        for source in ['PAO1_sources', 'PA14_sources', 'Ecoli_sources']:
            if filters[source][0] != 'None':
                sources.extend(filters[source])
        if len(sources) == 0:
            flash('Error: Please select sources from which to obtain interactions.')
            return redirect(request.url)

        filters['verification'] = [int(i) for i in filters['verification']]

        protein_file = None
        for root, dirs, files in os.walk(app.config['UPLOAD_FOLDER']):
            for file in files:
                filename = os.path.splitext(file)[0]
                if (filename == 'proteins'):
                    protein_file = os.path.join(app.config['UPLOAD_FOLDER'], file)

        input_proteins = []
        if protein_file is not None:
            with open(protein_file) as csvfile:
                delimiter = '\t'
                if protein_file.split('.')[-1] == 'csv':
                    delimiter = ','
                reader = csv.DictReader(csvfile, fieldnames=['ID'], delimiter = delimiter)
                for row in reader:
                    if row['ID'] != '':
                        input_proteins.append(row['ID'])

        session = Session()
        interactions = None
        if (len(filters['type']) == 1) & (filters['type'][0] == 'p-bs'):
            selected_tfbs = []
            if tfbs_sources[0] in sources:
                selected_tfbs.append(tfbs_sources[0])
            if tfbs_sources[1] in sources:
                selected_tfbs.append(tfbs_sources[1])
            if len(selected_tfbs) > 0:
                if len(input_proteins) == 0:
                    interactions = session.query(Interaction).filter(Interaction.strain.in_(filters['strain'])). \
                        join(Interaction.sources).filter(InteractionSource.data_source.in_(selected_tfbs),
                                                         InteractionSource.is_experimental.in_(
                                                             filters_all['verification'])).all()
                else:
                    interactions = session.query(Interaction).filter(Interaction.strain.in_(filters['strain'])). \
                        join(Interaction.sources, Interaction.interactors).\
                        filter(InteractionSource.data_source.in_(selected_tfbs),
                               InteractionSource.is_experimental.in_(filters_all['verification']),
                               Interactor.id.in_(input_proteins)).all()
        else:
            type = []
            if ('p-p' in filters['type']) | ('p-bs' in filters['type']):
                type.append('p-p')
            if 'p-m' in filters['type']:
                type.append('p-m')
                type.append('m-p')
            if len(input_proteins) == 0:
                interactions = session.query(Interaction).filter(
                    Interaction.strain.in_(filters['strain']), Interaction.type.in_(filters['type'])). \
                    join(Interaction.sources).filter(
                    InteractionSource.data_source.in_(sources),
                    InteractionSource.is_experimental.in_(filters['verification'])).all()
            else:
                interactions = session.query(Interaction).filter(
                    Interaction.strain.in_(filters['strain']), Interaction.type.in_(filters['type'])). \
                    join(Interaction.sources, Interaction.interactors).filter(
                    InteractionSource.data_source.in_(sources),
                    InteractionSource.is_experimental.in_(filters['verification']),
                    Interactor.id.in_(input_proteins)).all()

        file_writer = csv.DictWriter(open(
            (os.path.join(app.config['UPLOAD_FOLDER'], "output.csv")), mode='x', newline=''),
            fieldnames=psimi_fields, quoting=csv.QUOTE_NONE, delimiter='\t', quotechar=None)
        file_writer.writeheader()
        for interaction in interactions:
            if interaction is None: continue
            interactor_ids, alt_ids, aliases = [], [], []
            is_protein = []

            for interactor in interaction.interactors:
                if interactor.name is None:
                    aliases.append('-')
                else:
                    aliases.append(interactor.name)

                if interactor.type == 'p':
                    is_protein.append(1)
                    if interactor.uniprotkb == 'pc':
                        interactor_ids.append('uniprotkb:' + interactor.id)
                        alt_ids.append('-')
                    else:
                        interactor_ids.append('entrez_gene/locuslink:' + interactor.id)
                        alt_id = ''
                        if interactor.uniprotkb is not None:
                            alt_id += 'uniprotkb:' + interactor.uniprotkb + '|'
                        if interactor.ncbi_acc is not None:
                            alt_id += 'refseq:' + interactor.ncbi_acc + '|'
                        if len(alt_id) > 0:
                            alt_ids.append(alt_id[:-1])
                        else:
                            alt_ids.append('-')

                else:
                    is_protein.append(0)
                    id = None
                    alt_id = ''
                    if interactor.pubchem is not None:
                        id = 'pubchem:' + interactor.pubchem
                    if interactor.chebi is not None:
                        chebi = 'chebi:"CHEBI:' + interactor.chebi + '"'
                        if id is None:
                            id = chebi
                        else:
                            alt_id += chebi + '|'
                    if interactor.cas is not None:
                        cas = 'cas:' + interactor.cas
                        if id is None:
                            id = cas
                        else:
                            alt_id += cas + '|'
                    if interactor.kegg is not None:
                        kegg = 'kegg:' + interactor.kegg
                        if id is None:
                            id = kegg
                        else:
                            alt_id += kegg + '|'
                    if interactor.ecocyc is not None:
                        ecocyc = 'ecocyc:' + interactor.ecocyc
                        if id is None:
                            id = ecocyc
                        else:
                            alt_id += ecocyc + '|'

                    if len(alt_id) > 0:
                        alt_ids.append(alt_id[:-1])
                    else:
                        alt_ids.append('-')

            if len(interactor_ids) == 1:
                interactor_ids.append(interactor_ids[0])
                alt_ids.append(alt_ids[0])
                aliases.append(aliases[0])
                is_protein.append(is_protein[0])

            taxid_A, taxid_B, taxid = None, None, None
            if interaction.strain == 'PAO1':
                taxid = 'taxid:208964(pseae)'
            else:
                taxid = 'taxid:208963(pseab)'
            if is_protein[0]:
                taxid_A = taxid
            if is_protein[1]:
                taxid_B = taxid

            # 0064 ortholog interaction (interologs mapping)
            refs = {'detection': [], 'author': [], 'pmid': [], 'type': [], 'db': [], 'xrefs': [], 'confidence': [],
                    'annotations A': [], 'annotations B': []}
            author_temp = []
            for reference in interaction.references:
                refs['detection'].append(reference.detection_method)
                author_temp.append([reference.author_ln, reference.pub_date])
                refs['pmid'].append(reference.pmid)
                refs['type'].append(reference.interaction_type)
                refs['db'].append(reference.source_db)
                refs['confidence'].append(reference.confidence)
                refs['annotations A'].append(reference.interactor_a)
                refs['annotations B'].append(reference.interactor_b)

            for author_info in author_temp:
                author_full = ''
                if author_info[0] is not None:
                    author_full = author_info[0] + ' et al.'
                    if author_info[1] is not None:
                        author_full += ' (' + author_info[1] + ')'
                    refs['author'].append(author_full)
                else:
                    refs['author'].append(None)

            for xref in interaction.xrefs:
                if xref is not None:
                    refs['xrefs'].append(xref.data_source + ':' + xref.accession)
                else:
                    refs['xrefs'].append(xref)
            for ref in refs:
                refs[ref] = ['-' if field is None else field for field in refs[ref]]
                if len(refs[ref]) == 0:
                    refs[ref] = '-'
                    continue
                if refs[ref].count(refs[ref][0]) == len(refs[ref]):
                    refs[ref] = [refs[ref][0]]
                refs[ref] = '|'.join(refs[ref])

            file_writer.writerow({psimi_fields[0]: interactor_ids[0], psimi_fields[1]: interactor_ids[1],
                                  psimi_fields[2]: alt_ids[0], psimi_fields[3]: alt_ids[1],
                                  psimi_fields[4]: aliases[0], psimi_fields[5]: aliases[1],
                                  psimi_fields[6]: refs['detection'], psimi_fields[7]: refs['author'],
                                  psimi_fields[8]: refs['pmid'], psimi_fields[9]: taxid_A, psimi_fields[10]: taxid_B,
                                  psimi_fields[11]: refs['type'], psimi_fields[12]: refs['db'],
                                  psimi_fields[13]: refs['xrefs'], psimi_fields[14]: refs['confidence'],
                                  psimi_fields[15]: refs['annotations A'], psimi_fields[15]: refs['annotations B']})

        return render_template('query/download.html', filters=str(filters))
    return render_template('query/filter.html')

@bp.route('download')
def download():
    try:
        return send_file(os.path.join(app.config['UPLOAD_FOLDER'], "output.csv"),
                         attachment_filename="interactions.csv")
    except Exception as e:
        return str(e)


