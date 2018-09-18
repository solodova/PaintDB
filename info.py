from helpers import removeFiles
from flask import (
    Blueprint, render_template
)

bp = Blueprint('info', __name__, url_prefix='/info')
#app.register_blueprint(bp)

@bp.route('/data_sources')
def data_sources():
    removeFiles()
    return render_template('info/data_sources.html')