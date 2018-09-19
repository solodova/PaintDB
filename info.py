from helpers import removeFiles
from flask import (
    Blueprint, render_template
)
from flask import current_app as app

bp = Blueprint('info', __name__, url_prefix='/info')

@bp.route('/data_sources')
def data_sources():
    removeFiles(app)
    return render_template('info/data_sources.html')