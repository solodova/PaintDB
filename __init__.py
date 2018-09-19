# import os
# from flask import Flask
# from flask import render_template, redirect, url_for, request, flash
# from werkzeug.utils import secure_filename
# from sqlalchemy import create_engine
# from flask import current_app
# from sqlalchemy.orm import sessionmaker
# # from info import bp as bp_info
# # from query import bp as bp_query
# from helpers import removeFiles
#
#
#
# def create_app(test_config=None):
#     # create and configure the app
#     app = Flask(__name__, instance_relative_config=True)
#
#
#     app.config.from_mapping(
#         SECRET_KEY='dev',
#         UPLOAD_FOLDER=os.path.join('user_upload')
#     )
#     if test_config is None:
#         # load the instance config, if it exists, when not testing
#         app.config.from_pyfile('config.py', silent=True)
#     else:
#         # load the test config if passed in
#         app.config.from_mapping(test_config)
#
#     # ensure the instance folder exists
#     try:
#         os.makedirs(app.instance_path)
#     except OSError:
#         pass
#
#
#     # app.register_blueprint(bp_query)
#     # app.register_blueprint(bp_info)
#
#     # a simple page that says hello
#     @app.route('/')
#     def PaIntDB():
#         return redirect(url_for('home'))
#
#
#     @app.route('/home')
#     def home():
#         removeFiles()
#         return render_template('home.html')
#
#
#     return app