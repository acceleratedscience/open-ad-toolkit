import os
import json
import pandas
from flask import render_template, request, redirect, url_for
from ad4e_opentoolkit.helpers.output import output_table


def fetchRoutesDataViewer(data):
    def home():
        return render_template('/dataviewer/index.html', data=data)

    def submit():
        data_json = json.loads(request.data.decode('utf-8'))
        df = pandas.DataFrame(data_json)
        output_table(df, is_data=True)
        return 'Success'

    routes = {
        '/': {
            'func': home,
            'method': 'GET'
        },
        '/submit': {
            'func': submit,
            'method': 'POST'
        },
    }
    return routes
