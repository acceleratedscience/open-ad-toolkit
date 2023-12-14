"""performs search on a collection"""
import re
from copy import deepcopy
import readline
import numpy as np
from deepsearch.cps.client.components.elastic import ElasticDataCollectionSource
from deepsearch.cps.queries import DataQuery
from openad.helpers.output import output_table
from openad.helpers.output import output_text
from openad.helpers.output import output_error

# Importing our own plugins.
# This is temporary until every plugin is available as a public pypi package.
from openad.plugins.style_parser import style

_tableformat = "simple"

aggs = {
    "by_year": {
        "date_histogram": {
            "field": "description.publication_date",
            "calendar_interval": "year",
            "format": "yyyy",
            "min_doc_count": 0,
        }
    }
}


def search(inputs: dict, cmd_pointer):
    if cmd_pointer.notebook_mode is True:
        from tqdm.notebook import tqdm
        import IPython.display as display
    else:
        from tqdm import tqdm
    import pandas as pd

    search_query = ""
    val_index_key = "pubchem"
    page_size = 50
    val = "val"
    val_elastic_id = "default"
    highlight = {"fields": {"*": {}}}
    highlight["fragment_size"] = 0

    if "save_as" in inputs:
        highlight["pre_tags"] = [""]
        highlight["post_tags"] = [""]
    elif "return_as_data" in inputs:
        highlight["pre_tags"] = [""]
        highlight["post_tags"] = [""]
    elif cmd_pointer.notebook_mode is False:
        highlight["pre_tags"] = ["<green>"]
        highlight["post_tags"] = ["</green>"]
    else:
        highlight["pre_tags"] = ["<span style='font-weight: bold; background-color: #FFFF00'>"]
        highlight["post_tags"] = ["</span>"]

    edit_distance = 0
    search_query = inputs["search_text"]

    if "collection" in inputs:
        val_index_key = inputs["collection"]
    else:
        output_error("No collection_key suppled. ", cmd_pointer=cmd_pointer, return_val=False)
        return False
    if "elastic_id" in inputs:
        val_elastic_id = inputs["elastic_id"][val]
    if "page_size" in inputs:
        page_size = int(inputs["page_size"][val])
    if "edit_distance" in inputs:
        edit_distance = int(inputs["edit_distance"][val])
    else:
        edit_distance = 10
    if "display_first" in inputs:
        display_first = int(inputs["display_first"][val])
    else:
        display_first = 0

    api = cmd_pointer.login_settings["toolkits_api"][cmd_pointer.login_settings["toolkits"].index("DS4SD")]
    collections = api.elastic.list()
    collections.sort(key=lambda c: c.name.lower())

    elastic_list = [c.source.elastic_id for c in collections]

    index_list = [c.source.index_key for c in collections]
    index_name_list = [c.name for c in collections]

    result = [
        {
            "Domains": "/ ".join(c.metadata.domain),
            "Collection Name": c.name,
            "Collection key": c.source.index_key,
            "system_id": c.source.elastic_id,
        }
        for c in collections
    ]
    if val_elastic_id not in elastic_list:
        output_error("Invalid system_id, please choose from the following: ", cmd_pointer=cmd_pointer, return_val=False)
        collectives = pd.DataFrame(result)
        if cmd_pointer.notebook_mode is True:
            display.display(output_table(collectives, cmd_pointer=cmd_pointer))
        else:
            output_table(collectives, cmd_pointer=cmd_pointer)
        return False
    if val_index_key not in index_list and val_index_key not in index_name_list:
        output_error(
            "Invalid collection key or name, please choose from the following: ",
            cmd_pointer=cmd_pointer,
            return_val=False,
        )
        collectives = pd.DataFrame(result)

        if cmd_pointer.notebook_mode is True:
            display.display(output_table(collectives, cmd_pointer=cmd_pointer))
        else:
            output_table(collectives, cmd_pointer=cmd_pointer)
        return False

    if val_index_key in index_name_list:
        val_index_key = index_list[index_name_list.index(val_index_key)]

    data_collection = ElasticDataCollectionSource(elastic_id=val_elastic_id, index_key=val_index_key)

    # Prepare the data query
    source_list = []
    if edit_distance > 0:
        search_query = search_query + " ~" + str(edit_distance)

    is_docs = False

    if "show_data" in inputs:
        for i in inputs["show_data"]:
            if i.upper() == "DATA":
                source_list.extend(["subject", "attributes", "identifiers"])
            if i.upper() == "DOCS":
                source_list.extend(["description.title", "description.authors", "file-info.filename", "identifiers"])
                is_docs = True
    else:
        source_list = ["subject", "attributes", "identifiers", "file-info.filename"]

    if is_docs is False:
        edit_distance = 0

    if edit_distance > 0:
        query = DataQuery(
            search_query,  # The search query to be executed
            source=source_list,
            limit=page_size,  # The size of each request page
            highlight=highlight,
            coordinates=data_collection,  # The data collection to be queries
            aggregations=aggs,
        )
    else:
        query = DataQuery(
            search_query,  # The search query to be executed
            source=source_list,
            limit=page_size,  # The size of each request page
            coordinates=data_collection,  # The data collection to be queries
            aggregations=aggs,
        )

    # [Optional] Compute the number of total results matched. This can be used to monitor the pagination progress.
    count_query = deepcopy(query)
    count_query.paginated_task.parameters["limit"] = 0
    count_results = api.queries.run(count_query)

    expected_total = count_results.outputs["data_count"]
    expected_pages = (expected_total + page_size - 1) // page_size  # this is simply a ceiling formula

    if "estimate_only" in inputs:
        output_text("Expected Results Estimate: " + str(expected_total), cmd_pointer=cmd_pointer, return_val=False)
        return None
    else:
        output_text("\n Expected Results Estimate: " + str(expected_total), cmd_pointer=cmd_pointer, return_val=False)
        if expected_total > 100:
            if confirm_prompt("Your results may take some time to return, do you wish to proceed") is False:
                return None

    # Iterate through all results by fetching `page_size` results at the same time
    all_results = []
    all_aggs = {}
    try:
        cursor = api.queries.run_paginated_query(query)
    except Exception as e:  # pylint: disable=broad-exception-caught
        output_error("Error in calling deepsearch:" + str(e), cmd_pointer=cmd_pointer, return_val=False)
        return False

    for result_page in tqdm(cursor, total=expected_pages):
        all_results.extend(result_page.outputs["data_outputs"])
        for year in result_page.outputs["data_aggs"]["by_year"]["buckets"]:
            if year["key_as_string"] not in all_aggs:
                all_aggs[year["key_as_string"]] = 0
            all_aggs[year["key_as_string"]] = all_aggs[year["key_as_string"]] + int(year["doc_count"])

    results_table = []
    result = None

    if is_docs:
        df = pd.json_normalize(all_aggs)
        if cmd_pointer.notebook_mode is True:
            if len(df.columns) > 1:
                display.display(output_text("<h3>Distribution of Returned Documents by Year</h3>"))
                display.display(output_table(df, cmd_pointer=cmd_pointer))
        elif all_aggs != {}:
            if len(df.columns) > 1:
                output_text("\n<h1>Distribution of Returned Documents by Year</h1>")
                output_table(df, cmd_pointer=cmd_pointer)

    pd.set_option("display.max_colwidth", None)
    x = 0
    for row in all_results:
        while x < 20:
            x = x + 1
        result = {}
        if "description" in row["_source"]:
            if "title" in row["_source"]["description"]:
                result["Title"] = row["_source"]["description"]["title"]
            if "authors" in row["_source"]["description"]:
                result["Authors"] = ",".join([author["name"] for author in row["_source"]["description"]["authors"]])
            if "url_refs" in row["_source"]["description"]:
                result["URLs"] = " , ".join(row["_source"]["description"]["url_refs"])

        if edit_distance > 0:
            for field in row.get("highlight", {}).keys():
                for snippet in row["highlight"][field]:
                    result["Snippet"] = re.sub(" +", " ", snippet)

        if "attributes" in row["_source"]:
            for ref in row["_source"]["identifiers"]:
                if ref["type"] == "cid":
                    result["cid"] = ref["value"]
        for ref in row["_source"].get("identifiers", []):
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

        if "identifiers" in row["_source"]:
            for ref in row["_source"]["identifiers"]:
                if ref["type"] == "arxivid":
                    result["arxivid"] = f'https://arxiv.org/abs/{ref["value"]}'
                if ref["type"] == "doi":
                    result["doi"] = f'https://doi.org/{ref["value"]}'

        if edit_distance > 0:
            for field in row.get("highlight", {}).keys():
                for snippet in row["highlight"][field]:
                    result["Report"] = str(row["_source"]["file-info"]["filename"])
                    result["Field"] = field

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
        output_text("Search returned no result", cmd_pointer=cmd_pointer, return_val=False)
        return None
    if "save_as" in inputs:
        results_file = str(inputs["results_file"])
        if not results_file.endswith(".csv"):
            results_file = results_file + ".csv"
    if cmd_pointer.notebook_mode is True:
        df = pd.DataFrame(results_table)

        if "save_as" in inputs:
            df.to_csv(
                cmd_pointer.workspace_path(cmd_pointer.settings["workspace"].upper()) + "/" + results_file, index=False
            )
        df = df.replace(np.nan, "", regex=True)

        if "return_as_data" in inputs:
            return df
        else:
            df = df.style.format(hyperlinks="html")
            df = df.set_properties(**{"text-align": "left"})
            return df

    else:
        cmd_line_result = pd.DataFrame(results_table)

        if "save_as" in inputs:
            cmd_line_result.to_csv(
                cmd_pointer.workspace_path(cmd_pointer.settings["workspace"].upper()) + "/" + results_file, index=False
            )
        else:
            cmd_line_result.style.format(hyperlinks="html")
            if "Title" in cmd_line_result:
                cmd_line_result["Title"] = cmd_line_result["Title"].str.wrap(50, break_long_words=True)
            if "Authors" in cmd_line_result:
                cmd_line_result["Authors"] = cmd_line_result["Authors"].str.wrap(25, break_long_words=True)
            if "Snippet" in cmd_line_result:
                cmd_line_result["Snippet"] = cmd_line_result["Snippet"].apply(
                    lambda x: style(x)
                )  # pylint diable:unnecessary-lambda

                cmd_line_result["Snippet"] = cmd_line_result["Snippet"].str.wrap(70, break_long_words=True)
            if display_first > 0:
                cmd_line_result2 = cmd_line_result.truncate(after=display_first)
            else:
                cmd_line_result2 = cmd_line_result
            output_table(cmd_line_result2.replace(np.nan, "", regex=True), cmd_pointer, tablefmt=_tableformat)


def confirm_prompt(question: str) -> bool:
    reply = None
    while reply not in ("y", "n"):
        reply = input(f"{question} (y/n): ").casefold()
        readline.remove_history_item(readline.get_current_history_length() - 1)
    return reply == "y"
