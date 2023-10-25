

import os
import signal
import json
import pandas
from flask import render_template, request, redirect, url_for
from ad4e_opentoolkit.helpers.output import output_table


# TEMPORARY dummy cmd_pointer until we move memory functionality into its own class.
class TempDummyDmdPointer:
    notebook_mode = False
    memory = {}
    preserve_memory = {}
    api_mode = False


tempDummyDmdPointer = TempDummyDmdPointer()


def fetchRoutesDataViewer(data, notebook_mode):
    def home():
        return render_template('/dataviewer/index.html', data=data)

    def submit():
        data_json = json.loads(request.data.decode('utf-8'))
        df = pandas.DataFrame(data_json)
        output_table(df, tempDummyDmdPointer, is_data=True)  # return_value=False only for testing in Jupyter with tempDummyDmdPointer
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
