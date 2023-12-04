""" Displays a Domains Collections"""
import pandas as pd
from openad.helpers.output import output_table

_tableformat = "simple"


def display_collection(inputs: dict, cmd_pointer):
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
    if cmd_pointer.notebook_mode is True:
        return pd.DataFrame(results)
    else:
        collectives = pd.DataFrame(results)
        output_table(collectives, cmd_pointer, tablefmt=_tableformat)
