import functools, os
from flask import (
    Flask, Blueprint, flash, g, redirect, render_template, request, session, url_for
)
from werkzeug.utils import secure_filename
from flask import current_app as app

ALLOWED_EXTENSIONS = set(['txt', 'csv', 'tsv'])
bp = Blueprint('query', __name__, url_prefix='/query')
protein_file_name, metabolite_file_name = '', ''


def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS


@bp.route('upload', methods=['GET', 'POST'])
def upload():
    if request.method == 'POST':
        if ('proteinFile' not in request.files) & ('metaboliteFile' not in request.files):
            flash('Error: No files selected.')
            return redirect(request.url)
        protein_file, metabolite_file, metabolite_id = None, None, ''
        if ('metaboliteFile' in request.files):
            if 'metabolite_id' not in request.form:
                flash('Error: No ID type selected for metabolite interactor file.')
                return redirect(request.url)
            metabolite_file = request.files['metaboliteFile']
            metabolite_id = request.form['metabolite_id']
        if ('proteinFile' in request.files):
            protein_file = request.files['proteinFile']


        if (protein_file and not allowed_file(protein_file.filename)) or \
                (metabolite_file and not allowed_file(metabolite_file.filename)):
            if not allowed_file(protein_file.filename):
                flash('Error: File type for ' + protein_file.filename + ' is invalid')
            if not allowed_file(metabolite_file.filename):
                flash('Error: File type for ' + metabolite_file.filename + ' is invalid')
            redirect(request.url)


        if metabolite_file:
            if metabolite_file.filename == '':
                flash('Error: No file selected.')
                return redirect(request.url)
        if protein_file:
            if protein_file.filename == '':
                flash('Error: No file selected.')
                return redirect(request.url)

        if metabolite_file is not None:
            filename = secure_filename(metabolite_file)
            metabolite_file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
            flash(metabolite_file.filename + ' was successfully uploaded! Click below to move on to the next step.')
        if protein_file is not None:
            filename = secure_filename(protein_file)
            protein_file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
            flash(protein_file.filename + ' was successfully uploaded! Click below to move on to the next step.')

        # if user does not select file, browser also
        # submit an empty part without filename
        # if file.filename == '':
        #     flash('Error: No file selected.')
        #     return redirect(request.url)
        # if not allowed_file(file.filename):
        #     flash('Error: File type not supported. Please upload a txt, csv, or tsv file.')
        #     return redirect(request.url)
        # if file:
        #     filename = secure_filename(file.filename)
        #     file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
        #     protein_input = filename
        #     flash(file.filename + ' was successfully uploaded! Click below to move on to the next step.')
        #     return redirect(request.url)
    return render_template('query/upload.html')

@bp.route('upload/protein', methods=['GET', 'POST'])
def upload_protein():
    if request.method == 'POST':
        if 'file' not in request.files:
            flash('Error: No file selected.')
            return redirect(request.url)
        file = request.files['file']
        # if user does not select file, browser also
        # submit an empty part without filename
        if file.filename == '':
            flash('Error: No file selected.')
            return redirect(request.url)
        if not allowed_file(file.filename):
            flash('Error: File type not supported. Please upload a txt, csv, or tsv file.')
            return redirect(request.url)
        if file:
            filename = secure_filename(file.filename)
            file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
            protein_input = filename
            flash(file.filename + ' was successfully uploaded! Click below to move on to the next step.')
            return redirect(request.url)
    return render_template('query/uploadProtein.html')

@bp.route('upload/metabolite', methods=['GET', 'POST'])
def upload_metabolite():
    if request.method == 'POST':
        if 'metabolite_ids' not in request.form:
            flash('No ID type selected')
            return redirect(request.url)
        if 'file' not in request.files:
            flash('No selected file')
            return redirect(request.url)
        file = request.files['file']
        # if user does not select file, browser also
        # submit an empty part without filename
        if file.filename == '':
            flash('No selected file')
            return redirect(request.url)
        if not allowed_file(file.filename):
            flash('File type not supported. Please upload a txt, csv, or tsv file.')
            return redirect(request.url)
        if file:
            filename = secure_filename(file.filename)
            file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
            metabolite_input = filename
            flash(file.filename + ' was successfully uploaded! The ID type selected was: ' + \
                  request.form['metabolite_ids'] + '. Click below to move on to the next step.')
            return redirect(request.url)
    else:
        return render_template('query/uploadMetabolite.html')

@bp.route('filter', methods=['GET', 'POST'])
def filter():
    if request.method == 'POST':
        'hello'
    return render_template('query/filter.html')


