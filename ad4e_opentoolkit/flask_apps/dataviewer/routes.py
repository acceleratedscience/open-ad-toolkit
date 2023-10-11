import os
from flask import render_template


def fetchRoutesDataViewer(data):
    def home():
        return render_template('/dataviewer/index.html', data=data)

    routes = {
        '/': {
            'func': home,
            'method': 'GET'
        }
    }
    return routes
