import os
import signal
import json
import pandas
import threading
from flask import render_template, request
from openad.helpers.output import output_table
from openad.app.global_var_lib import GLOBAL_SETTINGS


def fetchRoutesDataViewer(data):
    def home():
        return render_template("/dataviewer/index.html", data=data)

    def submit():
        data_json = json.loads(request.data.decode("utf-8"))
        df = pandas.DataFrame(data_json)
        output_table(df)
        return "Success"

    def success():
        # Render success page before we kill the server, so we can still serve the CSS.
        html = render_template("/dataviewer/success.html", display=GLOBAL_SETTINGS["display"])

        # Kill server after 1 second so it has time to deliver the CSS etc.
        threading.Timer(1, lambda: os.kill(os.getpid(), signal.SIGINT)).start()

        return html

    routes = {
        "/": {"func": home, "method": "GET"},
        "/submit": {"func": submit, "method": "POST"},
        "/success": {"func": success, "method": "GET"},
    }
    return routes
