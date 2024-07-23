import json
import pandas
import IPython
from flask import Flask, jsonify, request


class GeneralApi:
    """
    All the general API endpoints.
    The API endpoints are called from gui_routes.py.
    """

    def __init__(self, cmd_pointer):
        self.cmd_pointer = cmd_pointer

    def landing(self):
        return "This is the OpenAD GUI API."

    def health(self):
        return ":)"

    def test(self):
        return "Hello World!"

    def exec_command(self):
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
            return "No result", 50
