import os
from flask import Flask
from flask import render_template, redirect, url_for

def create_app(test_config=None):

    # create and configure the app
    app = Flask(__name__, instance_relative_config=True)
    app.config.from_mapping(
        SECRET_KEY='dev',
        DATABASE=os.path.join(app.instance_path, 'flaskr.sqlite'),
    )

    if test_config is None:
        # load the instance config, if it exists, when not testing
        app.config.from_pyfile('config.py', silent=True)
    else:
        # load the test config if passed in
        app.config.from_mapping(test_config)

    # ensure the instance folder exists
    try:
        os.makedirs(app.instance_path)
    except OSError:
        pass

    # a simple page that says hello
    @app.route('/')
    def PaIntDB():
        return redirect(url_for('home'))

    @app.route('/home')
    def home():
        return render_template('home.html')

    from . import info
    app.register_blueprint(info.bp)

    from . import query
    app.register_blueprint(query.bp)

    # is blueprint needed for anything?
    # from . import query
    # app.register_blueprint(query.bp)

    return app