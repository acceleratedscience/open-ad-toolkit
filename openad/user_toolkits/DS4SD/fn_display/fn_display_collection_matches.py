# Example command:
# display collection matches for 'Ibuprofen'

import numpy as np
import pandas as pd
from deepsearch.cps.queries import DataQuery
from deepsearch.cps.client.components.queries import RunQueryError
from openad.helpers.output import output_error, output_table, output_success
from openad.helpers.output_msgs import msg
from openad.helpers.general import load_tk_module
from openad.app.global_var_lib import GLOBAL_SETTINGS


def display_collection_matches(inputs: dict, cmd_pointer):
    """
    Searches all collections for matches for a given string.

    Parameters
    ----------
    inputs:
        Parser inputs from pyparsing.
    cmd_pointer:
        Pointer to runtime.
    """

    # Load module from the toolkit folder.
    ds4sd_msg = load_tk_module(cmd_pointer, "DS4SD", "msgs", "ds4sd_msg")

    api = cmd_pointer.login_settings["toolkits_api"][cmd_pointer.login_settings["toolkits"].index("DS4SD")]
    if GLOBAL_SETTINGS["display"] == "notebook":
        from tqdm.notebook import tqdm
    else:
        from tqdm import tqdm
    try:
        collections = api.elastic.list()
        collections.sort(key=lambda c: c.name.lower())
        # raise Exception('This is a test error')
    except Exception as err:  # pylint: disable=broad-exception-caught
        output_error(ds4sd_msg("err_deepsearch", err), return_val=False)
        return False

    results = []
    for c in (pbar := tqdm(collections)):
        pbar.set_description(f"Querying {c.name}")

        # Search only on document collections
        if c.metadata.type != "Document":
            continue
        try:
            # Execute the query
            query = DataQuery(inputs["search_string"], source=[""], limit=0, coordinates=c.source)
            query_results = api.queries.run(query)
            if int(query_results.outputs["data_count"]) > 0:
                results.append(
                    {
                        "Domains": " / ".join(c.metadata.domain),
                        "Collection Name": c.name,
                        "Collection Key": c.source.index_key,
                        "Matches": query_results.outputs["data_count"],
                    }
                )
            # raise RunQueryError(task_id=1, message="This is a test error", error_type="err123", detail='aaa')
        except RunQueryError as err:
            output_error(ds4sd_msg("err_deepsearch", err), return_val=False)
            return False
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

    return output_table(pd.DataFrame(results))
