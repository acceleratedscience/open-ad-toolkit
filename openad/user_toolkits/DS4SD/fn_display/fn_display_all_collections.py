from openad.helpers.output import output_error, output_table, output_success
from openad.helpers.output_msgs import msg
from openad.app.global_var_lib import GLOBAL_SETTINGS

_tableformat = "simple"
import numpy as np


def display_all_collections(inputs: dict, cmd_pointer):
    """Displays all Collections
    inputs: parser inputs from pyparsing
    cmd_pointer: pointer to runtime"""
    import pandas as pd

    api = cmd_pointer.login_settings["toolkits_api"][cmd_pointer.login_settings["toolkits"].index("DS4SD")]
    try:
        collections = api.elastic.list()
    except Exception as err:  # pylint: disable=broad-exception-caught
        output_error(msg("err_deepsearch", err), return_val=False)
        return False

    collections.sort(key=lambda c: c.name.lower())
    results = [
        {
            "Domains": "/ ".join(c.metadata.domain),
            "Collection Name": c.name,
            "Collection key": c.source.index_key,
            "Type": c.metadata.type,
            "Num entries": c.documents,
            "Date": c.metadata.created.strftime("%Y-%m-%d"),
            "System": c.source.elastic_id,
        }
        for c in collections
    ]

    if "save_as" in inputs:
        results_file = str(inputs["results_file"])
        df = pd.DataFrame(results)
        if not results_file.endswith(".csv"):
            results_file = results_file + ".csv"
        df.to_csv(
            cmd_pointer.workspace_path(cmd_pointer.settings["workspace"].upper()) + "/" + results_file, index=False
        )
        df = df.replace(np.nan, "", regex=True)
        output_success(msg("success_file_saved"), return_val=False)

    if GLOBAL_SETTINGS["display"] == "notebook":
        return pd.DataFrame(results)
    else:
        collectives = pd.DataFrame(results)
        output_table(collectives, is_data=False, tablefmt=_tableformat)
