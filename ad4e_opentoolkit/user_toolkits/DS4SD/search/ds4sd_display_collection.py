import deepsearch as ds
import numerize
from deepsearch.cps.client.components.elastic import ElasticDataCollectionSource
from deepsearch.cps.queries import DataQuery
from deepsearch.cps.client.components.queries import RunQueryError
from  ad4e_opentoolkit.helpers.output import output_table as output_table
_tableformat = 'simple'


# login

# https://cps.foc-deepsearch.zurich.ibm.com/projects/1234567890abcdefghijklmnopqrstvwyz123456/library

def display_collection(inputs: dict, cmd_pointer):
    api = cmd_pointer.login_settings['toolkits_api'][cmd_pointer.login_settings['toolkits'].index('DS4SD') ]
    collections = api.elastic.list(domain=inputs['domain']['val'])
    collections.sort(key=lambda c: c.name.lower())

    results = [
        {
            "Name": c.name,
            "Type": c.metadata.type,
            "Num entries": c.documents,
            "Date": c.metadata.created.strftime("%Y-%m-%d"),
            "Coords": f"{c.source.elastic_id}/{c.source.index_key}",
        }
        for c in collections]
    if cmd_pointer.notebook_mode == True:

        import pandas as pd
        return pd.DataFrame(results)
    else:
        from tabulate import tabulate
        import pandas as pd
        collectives = pd.DataFrame(results)
        head = ["Name", "Type", "Num entries", "Date", "Coords"]
        #pdtabulate = tabulate(collectives, tablefmt=_tableformat, headers=head, showindex=False)
        #print('\n' + pdtabulate + '\n')
        output_table(collectives,cmd_pointer,tablefmt=_tableformat   )
           