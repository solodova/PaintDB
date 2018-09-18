import os
from flask import current_app

def removeFiles():
    for root, dirs, files in os.walk(current_app.config['UPLOAD_FOLDER']):
        for file in files:
            filename = os.path.splitext(file)[0]
            if (filename == 'proteins'):
                os.remove(os.path.join(current_app.config['UPLOAD_FOLDER'], file))