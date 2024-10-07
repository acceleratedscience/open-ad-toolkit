import json
import pandas as pd
from flask import request
from openad.helpers.files import open_file
from openad.helpers.output import output_table
from openad.app.global_var_lib import MEMORY
from openad.smols.smol_functions import (
    df_has_molecules,
    flatten_smol,
    create_molset_cache_file,
    assemble_cache_path,
    read_molset_from_cache,
)
from openad.smols.smol_transformers import dataframe2molset
from openad.gui.api.molecules_api import create_molset_response


class ResultApi:
    """
    All the API endpoints related to the result memory.
    The API endpoints are called from gui_routes.py.
    """

    def __init__(self, cmd_pointer):
        self.cmd_pointer = cmd_pointer

    def get_result(self):
        """
        Get data currently stored in result memory.
        """

        data = json.loads(request.data) if request.data else {}
        query = data["query"] if "query" in data else {}

        mem_data = MEMORY.get()

        # Nothing stored in memory.
        if mem_data is None:
            return {"type": "empty"}, 200

        # Memory has dataframe
        elif isinstance(mem_data, pd.DataFrame):
            # Dataframe has molecules -> load as molset.
            if df_has_molecules(mem_data):
                molset = dataframe2molset(mem_data)

                # Create cache working copy.
                cache_id = create_molset_cache_file(self.cmd_pointer, molset)

                # Read molset from cache.
                try:
                    molset = read_molset_from_cache(self.cmd_pointer, cache_id)
                except ValueError as err:
                    return f"get_result() -> {err}", 500

                return {"type": "molset", "data": create_molset_response(molset, query, cache_id)}, 200

            # Dataframe has no molecules -> load dataviewer.
            else:
                table = []
                for i, row in mem_data.iterrows():
                    mol = {}
                    for col in mem_data.columns:
                        mol[col] = None if pd.isna(row[col]) else row[col]
                    table.append(mol)

                # TO DO: Implement dataviewer here.
                return {"type": "data", "data": table}, 200

    def update_result_molset(self):
        """
        Save changes to a molset result stored in memory.
        """

        data = json.loads(request.data) if request.data else {}
        cache_id = data["cacheId"] if "cacheId" in data else ""

        if not cache_id:
            return f"update_result_molset() -> Unrecognized cache_id: {cache_id}", 500

        # Read data from cache.
        cache_path = assemble_cache_path(self.cmd_pointer, "molset", cache_id)
        molset, err_code = open_file(cache_path, return_err="code")
        if err_code:
            return err_code, 500

        # Flatten the mol dictionaries.
        molset = [flatten_smol(mol) for mol in molset]

        # Store the columns of the current result table,
        # so we can recreate them when overwriting the result.
        df = MEMORY.get()
        columns = df.columns.tolist()

        # Create new table.
        table = []
        for mol in molset:
            row = {}
            for i, col in enumerate(columns):
                row[col] = mol.get(col)
            table.append(row)

        # Write data back to memory as a dataframe.
        df = pd.DataFrame(table)
        MEMORY.store(df)

        return "ok", 200

    def update_result_data(self):
        """
        Placeholder for when we implement datavierwer for the /result page
        """
        return "ok", 200
