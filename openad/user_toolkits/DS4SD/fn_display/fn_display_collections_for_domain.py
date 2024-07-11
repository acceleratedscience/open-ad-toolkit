# Example command:
# display collections for domain 'Business Insights'

import pandas as pd
from openad.helpers.output import output_table


def display_collections_for_domain(inputs: dict, cmd_pointer):
    """
    Display available collections for a given domain.

    Parameters
    ----------
    inputs:
        Parser inputs from pyparsing.
    cmd_pointer:
        Pointer to runtime.
    """

    api = cmd_pointer.login_settings["toolkits_api"][cmd_pointer.login_settings["toolkits"].index("DS4SD")]
    collections = api.elastic.list(domain=inputs["domain"])
    collections.sort(key=lambda c: c.name.lower())
    results = [
        {
            "Name": c.name,
            "Type": c.metadata.type,
            "Num Entries": c.documents,
            "Date": c.metadata.created.strftime("%Y-%m-%d"),
            "Coords": f"{c.source.elastic_id}/{c.source.index_key}",
        }
        for c in collections
    ]

    return output_table(pd.DataFrame(results), is_data=False)
