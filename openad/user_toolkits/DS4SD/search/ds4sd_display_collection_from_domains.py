""" display all collections from a given Deep Sesrch Domain """
import numpy as np
import pandas as pd
from openad.helpers.output import output_table
from openad.helpers.output import output_error
from openad.helpers.output import output_text

_tableformat = "simple"


def display_all_collections(inputs: dict, cmd_pointer):
    """display all collections from a given Deep Sesrch Domain
    inputs: parser inputs from pyparsing
    cmd_pointer: pointer to runtime"""

    api = cmd_pointer.login_settings["toolkits_api"][cmd_pointer.login_settings["toolkits"].index("DS4SD")]
    try:
        collections = api.elastic.list()
    except Exception as e:  # pylint: disable=broad-exception-caught
        output_error("Error in calling deepsearch:" + str(e), cmd_pointer=cmd_pointer, return_val=False)
        return False

    collections.sort(key=lambda c: c.name.lower())
    from_list = []
    if isinstance(inputs["from_source"], dict) and inputs["from_source"]["from_list"] != None:
        try:
            from_list = inputs["from_source"]["from_list"]
        except Exception:  # pylint: disable=broad-exception-caught
            output_error(
                "unexpected pyparsing error. Please screenshot and report circumstance to OpenAD team",
                cmd_pointer=cmd_pointer,
                return_val=False,
            )
            output_error("Restart Notebook Kernel or application to proceed", cmd_pointer=cmd_pointer, return_val=False)
            return False

    elif "from_list" in inputs["from_source"][0]:
        try:
            from_list = inputs["from_source"][0]["from_list"]
        except Exception:  # pylint: disable=broad-exception-caught
            output_error(
                "unexpected pyparsing error. Please screenshot and report circumstance to OpenAD team",
                cmd_pointer=cmd_pointer,
                return_val=False,
            )
            output_error("Restart Notebook Kernel or application to proceed", cmd_pointer=cmd_pointer, return_val=False)
            return False

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

    filtered_list = []
    if len(from_list) > 0:
        for x in results:
            i = 0
            for y in from_list:
                if y.upper() in str(x["Domains"]).upper():
                    i = i + 1
            if i > 0:
                filtered_list.append(x)
        results = filtered_list

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

    if cmd_pointer.notebook_mode is True:
        return pd.DataFrame(results)
    else:
        collectives = pd.DataFrame(results)
        output_table(collectives, cmd_pointer, tablefmt=_tableformat)
