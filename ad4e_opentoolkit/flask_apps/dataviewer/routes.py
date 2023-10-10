import os
from flask import render_template


def fetchRoutesDataViewer(data, headers):
    def home():
        return render_template('/dataviewer/index.html', data=data, headers=headers)

    routes = {
        '/': {
            'func': home,
            'method': 'GET'
        }
    }
    return routes
