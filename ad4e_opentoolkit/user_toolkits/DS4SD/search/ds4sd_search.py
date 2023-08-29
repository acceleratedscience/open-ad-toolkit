import deepsearch as ds
from copy import deepcopy
from numerize import numerize
from deepsearch.cps.client.components.elastic import ElasticDataCollectionSource
from deepsearch.cps.queries import DataQuery
from deepsearch.cps.client.components.queries import RunQueryError
_tableformat = 'simple'


# login

# https://cps.foc-deepsearch.zurich.ibm.com/projects/1234567890abcdefghijklmnopqrstvwyz123456/library

def search(inputs: dict, toolkit_dir, cmd_pointer):
    if cmd_pointer.notebook_mode == True:
        from tqdm.notebook import tqdm
    else:
        from tqdm import tqdm

    search_query = ""
    # val_elastic_id     =   "default"
    val_index_key = "pubchem"
    page_size = 50
    val = "val"
    val_elastic_id = "default"

    search_query = inputs["search_text"]

    val_index_key = inputs["index_key"][val]

    if 'elastic_id' in inputs:
        val_elastic_id = inputs["elastic_id"][val]
    if 'page_size' in inputs:
        page_size = inputs["page_size"][val]

    data_collection = ElasticDataCollectionSource(elastic_id=val_elastic_id, index_key=val_index_key)

    api = cmd_pointer.login_settings['toolkits_api'][cmd_pointer.login_settings['toolkits'].index('DS4SD') ]
    # Prepare the data query
    source_list = []

    if "show_data" in inputs:
        for i in inputs["show_data"]:
            if i.upper() == "DATA":
                source_list.extend(["subject", "attributes", "identifiers"])
            if i.upper() == "DOCS":
                source_list.extend(["description.title", "description.authors", "identifiers"])
    else:
        source_list = ["subject", "attributes", "identifiers"]

    query = DataQuery(
        search_query,  # The search query to be executed
        # source=["subject", "attributes", "identifiers"], # Which fields of documents we want to fetch
        source=source_list,
        limit=page_size,  # The size of each request page
        coordinates=data_collection  # The data collection to be queries
    )

    # [Optional] Compute the number of total results matched. This can be used to monitor the pagination progress.
    count_query = deepcopy(query)
    count_query.paginated_task.parameters["limit"] = 0
    count_results = api.queries.run(count_query)

    expected_total = count_results.outputs["data_count"]
    expected_pages = (expected_total + page_size - 1) // page_size  # this is simply a ceiling formula

    if "estimate_only" in inputs:
        print("\n Expected Results Estimate: " + str(expected_total) + "\n")
        return None
    else:
        print("\n Expected Results Estimate: " + str(expected_total) + "\n")
        if expected_total > 100:
            if confirm_prompt("Your results may take some time to return, do you wish to proceed") == False:
                return None

    # Iterate through all results by fetching `page_size` results at the same time
    all_results = []
    cursor = api.queries.run_paginated_query(query)

    for result_page in tqdm(cursor, total=expected_pages):

        all_results.extend(result_page.outputs["data_outputs"])

    results_table = []
    result = None

    import pandas as pd
    # return pd.DataFrame.from_dict(pd.json_normalize(all_results), orient='columns')
    # return(pd.json_normalize(all_results))
    pd.set_option('display.max_colwidth', None)
    for row in all_results:

        result = {}

        if "description" in row["_source"]:

            if "title" in row["_source"]["description"]:
                result["Title"] = row["_source"]["description"]["title"]
            if "authors" in row["_source"]["description"]:

                # for author in row["_source"]["description"]["authors"]:
                result["Authors"] = ",".join([author["name"] for author in row["_source"]["description"]["authors"]])
                # result["Authors"]= atr
        if "attributes" in row["_source"]:
            for ref in row["_source"]["identifiers"]:
                if ref["type"] == "cid":
                    result["cid"] = ref["value"]
        for ref in row["_source"].get("identifiers", []):
            # if ref["type"] == "smiles":
            #     result["SMILES"] = ref["value"]
            result[ref["type"]] = ref["value"]
        if "subject" in row["_source"]:
            for ref in row["_source"]["subject"]["identifiers"]:
                if ref["type"] == "smiles":
                    result["SMILES"] = ref["value"]
                if ref["type"] == "echa_ec_number":
                    result["ec_number"] = ref["value"]
                if ref["type"] == "cas_number":
                    result["cas_number"] = ref["value"]
                if ref["type"] == "patentid":
                    result["Patent ID"] = ref["value"]
            for ref in row["_source"]["subject"]["names"]:
                if ref["type"] == "chemical_name":
                    result["chemical_name"] = ref["value"]

        if "attributes" in row["_source"]:
            for attribute in row["_source"]["attributes"]:
                for predicate in attribute["predicates"]:
                    value = predicate["value"]["name"]
                    if "nominal_value" in predicate:
                        value = predicate["nominal_value"]["value"]
                    elif "numerical_value" in predicate:
                        value = predicate["numerical_value"]["val"]
                    result[predicate["key"]["name"]] = value
        results_table.append(result)

    if result is None:
        print("Search returned no result")
        return None

    if cmd_pointer.notebook_mode == True:
        import pandas as pd
        df = pd.DataFrame(results_table)
        import numpy as np
        if 'save_as' in inputs:
            df.to_csv(cmd_pointer.workspace_path(cmd_pointer.settings['workspace'].upper()) + "/" + inputs['results_file'] + '.csv',)
        return (df. replace(np.nan, '', regex=True))

    else:
        from tabulate import tabulate
        import pandas as pd
        import numpy as np
        collectives = pd.DataFrame(results_table)

        if 'save_as' in inputs:

            collectives.to_csv(
                cmd_pointer.workspace_path(
                    cmd_pointer.settings['workspace'].upper()) +
                "/" +
                inputs['results_file'] +
                '.csv',
                index=False)
        else:
            head = ["Name", "Type", "Num entries", "Date", "Coords"]
            pdtabulate = tabulate(collectives.replace(np.nan, '', regex=True), tablefmt=_tableformat, headers="keys", showindex=False)
            print('\n' + pdtabulate + '\n')


def confirm_prompt(question: str) -> bool:
    import readline
    reply = None
    while reply not in ("y", "n"):
        reply = input(f"{question} (y/n): ").casefold()
        readline.remove_history_item(readline.get_current_history_length() - 1)
    return (reply == "y")
