import functools
import os
from flask import (
    Blueprint, flash, g, redirect, render_template, request, session, url_for, current_app as app
)

bp = Blueprint('info', __name__, url_prefix='/info')


def removeFiles():
    for root, dirs, files in os.walk(app.config['UPLOAD_FOLDER']):
        for file in files:
            filename = os.path.splitext(file)[0]
            if (filename == 'proteins') | (filename == 'metabolites'):
                os.remove(os.path.join(app.config['UPLOAD_FOLDER'], file))

@bp.route('/data_sources')
def data_sources():
    removeFiles()
    return render_template('info/data_sources.html')