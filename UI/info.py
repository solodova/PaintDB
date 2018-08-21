import functools

from flask import (
    Blueprint, flash, g, redirect, render_template, request, session, url_for
)

bp = Blueprint('info', __name__, url_prefix='/info')

@bp.route('/data_sources')
def data_sources():
    return render_template('info/data_sources.html')