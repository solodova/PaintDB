import functools, os, glob
from flask import (
    Flask, Blueprint, flash, g, redirect, render_template, request, url_for
)
from . import Session, removeFiles
from werkzeug.utils import secure_filename
from flask import current_app as app
from Schema import Interactor, Interaction, Protein, Metabolite, InteractionReference, InteractionXref, ProteinXref, \
InteractionSource, Base
import re, csv

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
    removeFiles()
    return render_template('query/upload.html')

@bp.route('uploaded', methods=['GET', 'POST'])
def uploaded():
    if request.method =='POST':
        removeFiles()
        return redirect(url_for('query.upload'))
    return render_template('query/uploaded.html')

@bp.route('filter', methods=['GET', 'POST'])
def filter():
    if request.method == 'POST':
        filters = {'strain': [], 'interaction_type': [], 'Ecoli_sources': [], 'PAO1_sources': [], 'PA14_sources': [],
                   'verification': []}
        sources_Ecoli = ['EcoCyc', 'RegulonDB(Ecoli)', 'IMEx(Ecoli)', 'BindingDB(Ecoli)', 'EBI-GOA-nonIntAct(Ecoli)',
                         'IntAct(Ecoli)', 'iRefIndex(Ecoli)', 'mentha(Ecoli)', 'MINT(Ecoli)', 'MPIDB(Ecoli)',
                         'UniProt(Ecoli)', 'DIP(Ecoli)', 'KEGG(Ecoli)']
        sources_PAO1 = ['Geoff', 'XLinkDB', 'Zhang', 'ADIPInteractomes(PAO1)', 'IMEx(PAO1)', 'IntAct(PAO1)',
                        'iRefIndex(PAO1)', 'mentha(PAO1)', 'MINT(PAO1)', 'Galan-Vasquez(PAO1)', 'KEGG(PAO1)']
        sources_PA14 = ['IMEx(PA14)', 'IntAct(PA14)', 'iRefIndex(PA14)', 'mentha(PA14)', 'MINT(PA14)', 'KEGG(PA14)',
                        'Galan-Vasquez(PA14)']
        filters_all ={'strain': ['PAO1', 'PA14'], 'interaction_type': ['p-p', 'p-m', 'p-bs'],
                      'Ecoli_sources': sources_Ecoli, 'PAO1_sources': sources_PAO1, 'PA14_sources': sources_PA14,
                      'verification': ['0', '1', '2']}
        sources = []
        tfbs_sources = ['RegulonDB(Ecoli)', 'Galan-Vasquez(PAO1)']


        for filter in filters:
            if filter in request.form:
                filters[filter] = request.form.getlist(filter)

        for filter in filters:
            if 'None' in filter:
                filters[filter]=['None']
            elif 'All' in filter:
                filters[filter]=filters_all[filter]

        if filters['PAO1_sources'][0] != 'None':
            sources.append(filters['PAO1_sources'])
        if filters['PA14_sources'][0] != 'None':
            sources.append(filters['PA14_sources'])
        if filters['Ecoli_sources'][0] != 'None':
            sources.append(filters['Ecoli_sources'])
        if len(sources) == (len(filters_all['PAO1_sources']) + len(filters_all['PA14_sources']) +
                            len(filters_all['Ecoli_sources'])):
            sources = ['All']
        if len(sources) == 0:
            flash('Error: Please select sources from which to obtain interactions.')
            return redirect(request.url)
        session = Session()

        interactions = None
        if len(filters['strain'] == 2):
            interactions = session.query(Interaction)
        elif 'PAO1' in filters['strain']:
            interactions = session.query(Interaction).filter_by(strain = 'PAO1')
        elif 'PA14' in filters['strain']:
            interactions = session.query(Interaction).filter_by(strain = 'PA14')

        type = []
        if (len(filters['interaction_type']) == 1) & (filters['interaction_type'][0] == 'bs'):
            selected_tfbs = []
            if tfbs_sources[0] in sources:
                selected_tfbs.append(tfbs_sources[0])
            if tfbs_sources[1] in sources:
                selected_tfbs.append(tfbs_sources[1])
            if len(selected_tfbs) > 0:
                interactions = interactions.join(InteractionSource).\
                    group_by(Interaction).filter(InteractionSource.data_source.in_(selected_tfbs))
        else:
            if ('p-p' in filters['interaction_type']) | ('p-bs' in filters['interaction_type']):
                type.append('p-p')
            if 'p-m' in filters['interaction_type']:
                type.append('p-m')
                type.append('m-p')
            if len(type) < 3:
                interactions = interactions.filter(Interaction.type.in_(type))

            interactions = interactions.join(Interaction.sources).group_by(Interaction)
            if sources[0] != 'All':
                interactions = interactions.filter(InteractionSource.data_source.in_(sources))

        if len(filters['verification']) < 3:
            interactions = interactions.filter(InteractionSource.is_experimental.in_(filters['verification']))

        protein_file = None
        for root, dirs, files in os.walk(app.config['UPLOAD_FOLDER']):
            for file in files:
                filename = os.path.splitext(file)[0]
                if (filename == 'proteins'):
                    protein_file = file

        if protein_file is not None:
            proteins = []
            with open(protein_file) as csvfile:
                delimiter = '\t'
                if protein_file.os.path.splitext(file)[1] == 'csv':
                    delimiter = ','
                reader = csv.DictReader(protein_file, fieldnames=['ID'], delimiter = delimiter)
                for row in reader:
                    if row['ID'] != '':
                        proteins.append(row['ID'])

            interactions = interactions.join('interactors').filter(Interaction.interactors.in_(proteins))


        interactions = interactions.join(Interaction.references)
        interactions = interactions.join(Interaction.xrefs)

        file_writer = csv.DictWriter(open('output.csv', mode='x', newline=''), fieldnames=psimi_fields)
        file_writer.writeheader()
        for interaction in interactions.all():
            if interaction is None: continue
            interactor_ids, alt_ids, aliases  = [], [], []
            is_protein = []

            for interactor in interaction.interactors:
                if interactor.name is None:
                    aliases.append('')
                else:
                    aliases.append(interactor.name)

                if interactor.type == 'p':
                    is_protein.append(1)
                    if interactor.uniprotkb == 'pc':
                        interactor_ids.append('uniprotkb:' + interactor.id)
                        alt_ids.append('')
                    else:
                        interactor_ids.append('gene/locus_link:' + interactor.id)
                        alt_id = ''
                        for xref in interactor.xrefs:
                            alt_id += xref.source + ':' + xref.accession + '|'
                        if len(alt_id) != 0:
                            alt_id = alt_id[:-1]
                        alt_ids.append(alt_id)

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

                    if len(alt_id) != 0:
                        alt_id = alt_id[:-1]

                    alt_ids.append(alt_id)

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
            refs = {'detection': [], 'author': [], 'pmid': [], 'type': [], 'db': [], 'xrefs': [],
                    'confidence': [], 'annotations A': [], 'annotations B': []}
            #sorted(list_with_none, key=lambda k: (k[col] is not None, k[col] != "", k[col]), reverse=True)
            for reference in interaction.references:
                refs['detection'].append(reference.detection)
                refs['author'].append([reference.author_ln, reference.pub_date])
                refs['type'].append(reference.interaction_type)
                refs['db'].append(reference.source_db)
                refs['confidence'].append(reference.confidence)
                refs['annotations A'].append(reference.interactor_a)
                refs['annotations B'].append(reference.interactor_b)

            for author_info in refs['author']:
                if author_info[0] is not None:
                    author_info[0] += 'et al.'
                    if author_info[1] is not None:
                        author_info[0] += ' (' + author_info[1] + ')'
                refs['author'][author_info] = author_info[0]

            for xref in interaction.xrefs:
                if xref is None: continue
                refs['xrefs'].append(xref.data_source + ':' + xref.accession)
            for ref in refs:
                refs[ref] = ['-' if field is None else field for field in refs[ref]]
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

        return render_template('query/results.html', filters=str(filters))
    return render_template('query/filter.html')


