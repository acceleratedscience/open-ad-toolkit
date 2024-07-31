import os
import signal
import json
import pandas
import threading
from flask import render_template, request
from openad.helpers.output import output_table
from openad.app.global_var_lib import GLOBAL_SETTINGS
import openad.helpers.jupyterlab_settings as jl_settings


JL_PROXY = False
IS_STATIC = ""
try:
    jl = jl_settings.get_jupyter_lab_config()

    if (
        jl["ServerApp"]["allow_remote_access"] is True
        and "127.0.0.1" in jl["ServerProxy"]["host_allowlist"]
        and GLOBAL_SETTINGS["display"] == "notebook"
    ):
        JL_PROXY = True
        IS_STATIC = "/static"
except Exception as e:
    JL_PROXY = False
    IS_STATIC = ""


def fetchRoutesDataViewer(data):
    def home():
        return render_template(f"/dataviewer/{IS_STATIC}/index.html", data=data)

    def submit():
        data_json = json.loads(request.data.decode("utf-8"))
        df = pandas.DataFrame(data_json)
        output_table(df)
        return "Success"

    def success():
        # Render success page before we kill the server, so we can still serve the CSS.
        html = render_template(f"/dataviewer/{IS_STATIC}/success.html", display=GLOBAL_SETTINGS["display"])
        # Kill server after 1 second so it has time to deliver the CSS etc.
        threading.Timer(1, lambda: os.kill(os.getpid(), signal.SIGINT)).start()

        return html

    routes = {
        "/": {"func": home, "method": "GET"},
        "/submit": {"func": submit, "method": "POST"},
        "/success": {"func": success, "method": "GET"},
    }
    return routes
