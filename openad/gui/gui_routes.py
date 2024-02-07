"""
This file contains all the API endpoints consumed by the GUI.
"""

import json
from flask import request
from openad.core.lang_file_system import get_workspace_files as _get_workspace_files


def fetchRoutes(cmd_pointer):
    api_v1 = "/api/v1"

    def test():
        return "Hello World!"

    # Fetch your active workspace's content as a JSON object.
    def get_workspace_files():
        # data = request.data.decode("utf-8")
        data = json.loads(request.data) if request.data else {}
        path = data["path"] if "path" in data else ""
        return _get_workspace_files(cmd_pointer, path)

    routes = {
        f"{api_v1}/get-workspace-files": {"func": get_workspace_files, "method": "POST"},
        f"{api_v1}/test": {"func": test, "method": "GET"},
        # "/": {"func": home, "method": "GET"},
        # "/submit": {"func": submit, "method": "POST"},
        # "/success": {"func": success, "method": "GET"},
    }
    return routes
