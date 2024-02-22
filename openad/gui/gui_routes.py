"""
This file contains all the API endpoints consumed by the GUI.

They can be accessed on: http://localhost:5000/api/v1/<endpoint>
"""

import json
import pandas
import IPython
from flask import request
from openad.core.lang_file_system import get_workspace_files as _get_workspace_files
from openad.core.lang_file_system import get_file as _get_file


def fetchRoutes(cmd_pointer):
    api_v1 = "/api/v1"

    def test():
        return "Hello World!"

    # Fetch the name of the active workspace.
    def get_workspace_name():
        return cmd_pointer.settings["workspace"].upper()

    # Fetch your active workspace's content as a JSON object.
    def get_workspace_files():
        # data = request.data.decode("utf-8")
        data = json.loads(request.data) if request.data else {}
        path = data["path"] if "path" in data else ""
        return _get_workspace_files(cmd_pointer, path)

    # Fetch a file's content as a JSON object.
    def get_file():
        data = json.loads(request.data) if request.data else {}
        path = data["path"] if "path" in data else ""
        return _get_file(cmd_pointer, path)

    def exec_command():
        from openad.app.main import api_remote

        data = json.loads(request.data) if request.data else {}
        command = data["command"] if "command" in data else ""
        # print("Parsing command:\n", command)
        response = api_remote(command)
        # print("Response:\n", response)

        if hasattr(pandas.io.formats, "style") and isinstance(response, pandas.io.formats.style.Styler):
            # print("Response is pandas Styler object")
            response = response.data

        if isinstance(response, IPython.core.display.Markdown):
            # print("Response is IPython Markdown object")
            from IPython.core.formatters import format_display_data

            formatted, metadata = format_display_data(response)
            response = formatted["text/markdown"]

        if isinstance(response, pandas.core.frame.DataFrame):
            # print("Response is pandas DataFrame")
            response = response.to_csv()

        if response:
            return response, 200
        else:
            return "No result", 500

    routes = {
        f"{api_v1}/test": {"func": test, "method": "GET"},
        f"{api_v1}/get-workspace-name": {"func": get_workspace_name, "method": "GET"},
        f"{api_v1}/get-workspace-files": {"func": get_workspace_files, "method": "POST"},
        f"{api_v1}/get-file": {"func": get_file, "method": "POST"},
        f"{api_v1}/exec-command": {"func": exec_command, "method": "POST"},
        # "/": {"func": home, "method": "GET"},
        # "/submit": {"func": submit, "method": "POST"},
        # "/success": {"func": success, "method": "GET"},
    }
    return routes
