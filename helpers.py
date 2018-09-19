import os


def removeFiles(app):
    for root, dirs, files in os.walk(app.config['UPLOAD_FOLDER']):
        for file in files:
            filename = os.path.splitext(file)[0]
            if (filename == 'proteins') | (filename == 'output'):
                os.remove(os.path.join(app.config['UPLOAD_FOLDER'], file))