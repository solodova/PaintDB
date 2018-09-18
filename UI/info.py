import functools
import os
from app import removeFiles
from flask import (
    Blueprint, flash, g, redirect, render_template, request, session, url_for, current_app as app
)
from flask import current_app as app

bp = Blueprint('info', __name__, url_prefix='/info')
app.register_blueprint(bp)

@bp.route('/data_sources')
def data_sources():
    removeFiles()
    return render_template('info/data_sources.html')