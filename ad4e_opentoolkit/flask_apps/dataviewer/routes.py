

import os
import signal
import json
import pandas
from flask import render_template, request, redirect, url_for
from ad4e_opentoolkit.helpers.output import output_table



def fetchRoutesDataViewer(data, cmd_pointer):
    notebook_mode = cmd_pointer.notebook_mode

    def home():
        return render_template('/dataviewer/index.html', data=data)

    def submit():
        data_json = json.loads(request.data.decode('utf-8'))
        df = pandas.DataFrame(data_json)
        output_table(df, cmd_pointer, is_data=True)
        return 'Success'

    def success():
        # Kill server
        os.kill(os.getpid(), signal.SIGINT)

        return render_template('/dataviewer/success.html', notebook_mode=notebook_mode)

    routes = {
        '/': {
            'func': home,
            'method': 'GET'
        },
        '/submit': {
            'func': submit,
            'method': 'POST'
        },
        '/success': {
            'func': success,
            'method': 'GET'
        },
    }
    return routes
