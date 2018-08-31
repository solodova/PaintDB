import functools, os, glob
from flask import (
    Flask, Blueprint, flash, g, redirect, render_template, request, url_for
)
from . import Session, removeFiles
from werkzeug.utils import secure_filename
from flask import current_app as app
from Schema1 import Interactor, Interaction, Protein, Metabolite, InteractionReference, InteractionXref, InteractorXref, \
InteractionSource
import re
ALLOWED_EXTENSIONS = set(['txt', 'csv', 'tsv'])
bp = Blueprint('query', __name__, url_prefix='/query')
import sqlalchemy

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

@bp.route('upload', methods=['GET', 'POST'])
def upload():
    if request.method == 'POST':
        # if no files were added before submit, flash error
        if ('proteinFile' not in request.files) & ('metaboliteFile' not in request.files):
            flash('Error: No files selected.')
            return redirect(request.url)

        protein_file, metabolite_file, metabolite_id = None, None, ''
        file_error = ''
        # if a metabolite file was included but no id was selected, flash error
        if 'metaboliteFile' in request.files:
            if 'metabolite_id' not in request.form:
                flash('Error: No ID type selected for metabolite interactor file.')
                return redirect(request.url)
            metabolite_file = request.files['metaboliteFile']
            metabolite_id = request.form['metabolite_id']
            if metabolite_file.filename == '':
                flash('Error: No files selected.')
                return redirect(request.url)
            if not allowed_file(metabolite_file.filename):
                file_error += metabolite_file.filename + ' '

        if 'proteinFile' in request.files:
            protein_file = request.files['proteinFile']
            if protein_file.filename == '':
                flash('Error: No files selected.')
                return redirect(request.url)
            if not allowed_file(protein_file.filename):
                file_error += protein_file.filename + ' '

        if len(file_error) != 0:
            flash('Error: Incompatible file type for: ' + file_error)
            return redirect(request.url)

        if metabolite_file is not None:
            mfile = secure_filename(('metabolites.' + metabolite_file.filename.rsplit('.', 1)[1].lower()))
            metabolite_file.save(os.path.join(app.config['UPLOAD_FOLDER'], mfile))
            flash(metabolite_file.filename + ' was successfully uploaded! Click below to move on to the next step.')
        if protein_file is not None:
            #pfile = secure_filename(('proteins.' + protein_file.filename.rsplit('.', 1)[1].lower()))
            pfile = 'proteins.' + protein_file.filename.rsplit('.', 1)[1].lower()
            parts = re.split(r'/', pfile)
            protein_file.save(os.path.join(app.config['UPLOAD_FOLDER'], *parts))
            #protein_file.save('UI\\user_upload\\' + pfile)
            flash(protein_file.filename + ' was successfully uploaded! Click below to move on to the next step.')
        return redirect(url_for('query.uploaded'))
    removeFiles()
    return render_template('query/upload.html')

@bp.route('uploaded', methods=['GET', 'POST'])
def uploaded():
    if request.method =='POST':
        for root, dirs, files in os.walk(app.config['UPLOAD_FOLDER']):
            for file in files:
                filename = os.path.splitext(file)[0]
                if (filename == 'proteins') | (filename == 'metabolites'):
                    os.remove(os.path.join(app.config['UPLOAD_FOLDER'], file))
        return redirect(url_for('query.upload'))
    return render_template('query/uploaded.html')

@bp.route('filter', methods=['GET', 'POST'])
def filter():
    if request.method == 'POST':
        filters = {'strain': [], 'interaction_type': [], 'ortholog_mapping': [], 'Ecoli_sources': [],
                   'PAO1_sources': [], 'PA14_sources': [], 'verification': []}
        sources_Ecoli = ['EcoCyc', 'RegulonDB', 'IMEx(Ecoli)', 'BindingDB(Ecoli)', 'EBI-GOA-nonIntAct(Ecoli)',
                         'IntAct(Ecoli)',
                         'iRefIndex(Ecoli)', 'mentha(Ecoli)', 'MINT(Ecoli)', 'MPIDB(Ecoli)', 'UniProt(Ecoli)',
                         'DIP(Ecoli)', 'KEGG(Ecoli)']
        sources_PAO1 = ['Geoff', 'XLinkDB', 'Zhang', 'ADIPInteractomes(PAO1)', 'IMEx(PAO1)', 'IntAct(PAO1)',
                        'iRefIndex(PAO1)', 'mentha(PAO1)', 'MINT(PAO1)', 'Galan-Vasquez(PAO1)', 'KEGG(PAO1)']
        sources_PA14 = ['IMEx(PA14)', 'IntAct(PA14)', 'iRefIndex(PA14)', 'mentha(PA14)', 'MINT(PA14)', 'KEGG(PA14)',
                        'Galan-Vasquez(PA14)']
        filters_all ={'strain': ['PAO1', 'PA14'], 'interaction_type': ['p-p', 'p-m', 'p-bs'],
                      'ortholog_mapping': ['PAO1/PA14', 'PAO1/Ecoli', 'PA14/Ecoli'], 'Ecoli_sources': sources_Ecoli,
                      'PAO1_sources': sources_PAO1, 'PA14_sources': sources_PA14, 'verification': ['0, 1, 2']}
        sources = []
        session = Session()

        for filter in filters:
            if filter in request.form:
                filters[filter] = request.form.getlist(filter)

        for filter in filters:
            if 'None' in filter:
                filters[filter]=['None']
            elif 'All' in filter:
                filters[filter]=['All']
            elif len(filters[filter]) == len(filters_all[filter]):
                filters[filter] = ['All']

        if ('PAO1' in filters['strain']) | ('PAO1/PA14' in filters['ortholog_mapping']):
            sources.append(filters['PAO1_sources'])
        if ('PA14' in filters['strain']) | ('PAO1/PA14' in filters['ortholog_mapping']):
            sources.append(filters['PA14_sources'])
        if ('PAO1/Ecoli' in filters['ortholog_mapping']) | ('PA14/Ecoli' in filters['ortholog_mapping']):
            sources.append(filters['Ecoli_sources'])
        session = Session()

        interactions = None
        if (('PAO1' in filters['strain']) and ('PA14' in filters['strain'])) | ('All' in filters['strain']):
            interactions = session.query(Interaction)
        elif 'PAO1' in filters['strain']:
            interactions = session.query(Interaction).filter_by(strain = 'PAO1')
        elif 'PA14' in filters['strain']:
            interactions = session.query(Interaction).filter_by(strain = 'PAO1')

        interactions = session.query(Interaction).join(Interaction.sources)

        if filters['strain'][0] != 'All':
            interactions = interactions.filter(Interaction.strain.in_(filters['strain']))
        if filters['interaction_type'] != 'All':
            types = []
            if 'p-p' in filters['interaction_type']:
                types.append('p-p'), types.append('p-pc'), types.append('pc-p')
            if 'p-m' in filters['interaction_type']:
                types.append('p-m'), types.append('m-p')
            if len(types) != 0:
                interactions = interactions.filter(Interaction.type.in_(filters['interaction_type']))
        if sources



        #for strain in filters['strain']:
            #session.query(Interaction).filter(Interaction.strain == strain)

        return render_template('query/results.html', filters=str(filters))
    return render_template('query/filter.html')


