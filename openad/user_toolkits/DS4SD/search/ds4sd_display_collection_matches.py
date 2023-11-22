""" performs a search of collection for a specified string"""
import numpy as np
import pandas as pd
from deepsearch.cps.queries import DataQuery
from deepsearch.cps.client.components.queries import RunQueryError
from openad.helpers.output import output_table
from openad.helpers.output import output_error
from openad.helpers.output import output_text

_tableformat = "simple"


def display_collection_matches(inputs: dict, cmd_pointer):
    """Searches all collections for matches for a given Search String
    inputs: parser inputs from pyparsing
    cmd_pointer: pointer to runtime
    """

    api = cmd_pointer.login_settings["toolkits_api"][cmd_pointer.login_settings["toolkits"].index("DS4SD")]
    if cmd_pointer.notebook_mode is True:
        from tqdm.notebook import tqdm
    else:
        from tqdm import tqdm
    try:
        collections = api.elastic.list()
        collections.sort(key=lambda c: c.name.lower())
    except Exception as e:  # pylint: disable=broad-exception-caught
        output_error("Error in calling deepsearch:" + str(e), cmd_pointer=cmd_pointer, return_val=False)
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
                        "Domains": "/ ".join(c.metadata.domain),
                        "Collection Name": c.name,
                        "Collection Key": c.source.index_key,
                        "matches": query_results.outputs["data_count"],
                    }
                )
        except RunQueryError as err:
            output_error("Error in callling deepsearch:" + str(err), cmd_pointer=cmd_pointer, return_val=False)
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
        output_text(
            "\n <success>File successfully saved to workspace.</success>", cmd_pointer=cmd_pointer, return_val=False
        )

    collectives = pd.DataFrame(results)
    return output_table(collectives, cmd_pointer, tablefmt=_tableformat)
