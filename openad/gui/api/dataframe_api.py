import json
import pandas as pd
from flask import jsonify, request


class DataframeApi:
    """
    All the API endpoints related to dataframes referred to in Jupyter magic commands.
    The API endpoints are called from gui_routes.py.
    """

    def __init__(self, cmd_pointer):
        self.cmd_pointer = cmd_pointer

    def get_dataframe(self, df_name):
        """
        fetch the data from a dataframe referred to in a magic command.
        """

        if self.cmd_pointer.api_variables:
            df = self.cmd_pointer.api_variables.get(df_name)
            if df is not None and not df.empty:
                json_data = df.to_json(orient="records")
                if json_data is not None:
                    return jsonify(json_data), 200

        return {"type": "empty"}, 200

    def update_dataframe(self):
        """
        Save changes to a datraframe referred to in a magic command.
        """

        return "ok", 200
