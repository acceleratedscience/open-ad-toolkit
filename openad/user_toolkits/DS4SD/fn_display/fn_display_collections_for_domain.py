""" Displays a Domains Collections"""
import pandas as pd
from openad.helpers.output import output_table
from openad.app.global_var_lib import GLOBAL_SETTINGS

_tableformat = "simple"


def display_collections_for_domain(inputs: dict, cmd_pointer):
    """Displays a Domains Collections"""
    api = cmd_pointer.login_settings["toolkits_api"][cmd_pointer.login_settings["toolkits"].index("DS4SD")]
    collections = api.elastic.list(domain=inputs["domain"])
    collections.sort(key=lambda c: c.name.lower())
    results = [
        {
            "Name": c.name,
            "Type": c.metadata.type,
            "Num entries": c.documents,
            "Date": c.metadata.created.strftime("%Y-%m-%d"),
            "Coords": f"{c.source.elastic_id}/{c.source.index_key}",
        }
        for c in collections
    ]
    if GLOBAL_SETTINGS["display"] == "notebook":
        return pd.DataFrame(results)
    else:
        collectives = pd.DataFrame(results)
        output_table(collectives, tablefmt=_tableformat)
