# Example command:
# display all collections

import numpy as np
from openad.helpers.output import output_error, output_table, output_success
from openad.helpers.output_msgs import msg
from openad.helpers.general import load_tk_module


def display_all_collections(inputs: dict, cmd_pointer):
    """
    Display all collections.

    Parameters
    ----------
    inputs:
        Parser inputs from pyparsing.
    cmd_pointer:
        Pointer to runtime.
    """

    import pandas as pd

    # Load module from the toolkit folder.
    ds4sd_msg = load_tk_module(cmd_pointer, "DS4SD", "msgs", "ds4sd_msg")

    api = cmd_pointer.login_settings["toolkits_api"][cmd_pointer.login_settings["toolkits"].index("DS4SD")]
    try:
        collections = api.elastic.list()
        # raise Exception('This is a test error')
    except Exception as err:  # pylint: disable=broad-exception-caught
        output_error(ds4sd_msg("err_deepsearch", err), return_val=False)
        return False

    collections.sort(key=lambda c: c.name.lower())
    results = [
        {
            "Domains": " / ".join(c.metadata.domain),
            "Collection Name": c.name,
            "Collection Key": c.source.index_key,
            "Type": c.metadata.type,
            "Num Entries": c.documents,
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
        output_success(msg("success_file_saved", results_file), return_val=False, pad_top=1, pad_btm=0)

    return output_table(pd.DataFrame(results), is_data=False)
