# PaintDB

Note to new users:

- PaIntDB was built using sqlalchemy, and the web app was created using Flask. Sqlalchemy and Flask are python packages.

- DB_build.py contains code to create the database and to get the database statistics. 
This is the main executable for the project. If you wish to run specific parts (eg. just the statistics),
comment out the instructions you are not interested in. NOTE: If you run instructions that create files, 
make sure to delete/rename existing files (eg. for DB_statistics operations, delete/rename files in 
DB_statistics/PaIntDB_stats, to rebuild the whole database, delete/rename the current database file you have saved).

- The route to the PaIntDB.db file may need to be changed to match the route on your computer. The two places in the project where this route is specified are: DB_build.py (line 11) and query.py (line 11). If rebuilding the database or getting statistics, fix the route in DB_build.py, to use the web app, fix the route in query.py. NOTE: PaIntDB is an SQLite database, so use the approproate configuration. Refer to the sqlalchemy documentation on engine configuration (http://docs.sqlalchemy.org/en/latest/core/engines.html, relevant section is copied below) for the OS you are operating on.
    - relevant section for engine configuration:
        - Unix/Mac - 4 initial slashes in total
            - engine = create_engine('sqlite:////absolute/path/to/foo.db')
        - Windows
            - engine = create_engine('sqlite:///C:\\path\\to\\foo.db')
        - Windows alternative using raw string
            - engine = create_engine(r'sqlite:///C:\path\to\foo.db')

- to use the web app:
    - from terminal, cd to the project directory.
    - when running ls command, APP.py should show up in the project folder
    - run the following commands (NOTE: on Windows, use 'set' instead of 'export'):
          export FLASK_APP=APP.py
          export FLASK_ENV=development
          flask run
          
    - refer to http://flask.pocoo.org/docs/1.0/quickstart/ for additional documentation on running the web app.
    
