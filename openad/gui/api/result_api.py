import json
import pandas as pd
from flask import request
from openad.helpers.files import open_file
from openad.helpers.output import output_table
from openad.app.global_var_lib import MEMORY
from openad.molecules.mol_functions import (
    df_has_molecules,
    molformat_v2_to_v1,
    create_molset_cache_file,
    assemble_cache_path,
    read_molset_from_cache,
)
from openad.molecules.mol_transformers import dataframe2molset
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

    def update_molset_result(self):
        """
        Save changes to a molset result stored in memory.
        """
        data = json.loads(request.data) if request.data else {}
        cache_id = data["cacheId"] if "cacheId" in data else ""

        if not cache_id:
            return f"update_result() -> Unrecognized cache_id: {cache_id}", 500

        # Read data from cache.
        cache_path = assemble_cache_path(self.cmd_pointer, "molset", cache_id)
        molset, err_code = open_file(cache_path, return_err="code")
        if err_code:
            return err_code, 500

        # Flatten the mol dictionaries.
        molset = [molformat_v2_to_v1(mol) for mol in molset]
        props = set()
        for mol in molset:
            props.update(mol["properties"])

        # Store the columns of the current result table,
        # so we can recreate them when overwriting the result.
        df = MEMORY.get()
        columns = df.columns.tolist()
        columns_lower = [col.lower() for col in columns]  # Lets us match case-insensitive

        # Create new table
        table = []
        for mol in molset:
            row = {}
            for i, col_lower in enumerate(columns_lower):
                col = columns[i]
                row[col] = mol["properties"][col_lower]
            table.append(row)

        # Write data back to memory as a dataframe.
        df = pd.DataFrame(table)
        MEMORY.store(df)

        # Print result
        output_table(df)

        return "ok", 200

    def update_data_result(self):
        """
        Placeholder for when we implement datavierwer for the /result page
        """
        return "ok", 200
