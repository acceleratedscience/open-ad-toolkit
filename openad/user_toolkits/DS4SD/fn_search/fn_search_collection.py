# Example commands:
# search collection 'arxiv-abstract' for 'ide(\"power conversion efficiency\" OR PCE) AND organ*' using ( edit_distance=20 system_id=default) show (docs)
# search collection 'pubchem' for 'Ibuprofen' show (data)
# search collection 'patent-uspto' for '\"smiles#ccc(coc(=o)cs)(c(=o)c(=o)cs)c(=o)c(=o)cs\"' show (data)

import re
import base64
import json
import urllib.parse
import os
from copy import deepcopy
import readline
import numpy as np
from deepsearch.cps.client.components.elastic import ElasticDataCollectionSource, ElasticProjectDataCollectionSource
from typing import TYPE_CHECKING, Any, Dict, Literal, Optional, Union
from deepsearch.cps.queries import DataQuery
from openad.app.global_var_lib import GLOBAL_SETTINGS
from openad.plugins.style_parser import style
from openad.helpers.output import output_text, output_table, output_error
from openad.helpers.credentials import load_credentials
from openad.helpers.general import load_tk_module

DEFAULT_URL = "https://sds.app.accelerate.science/"


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


def search_collection(inputs: dict, cmd_pointer):
    """
    Perform search on a collection.

    Parameters
    ----------
    inputs:
        Parser inputs from pyparsing.
    cmd_pointer:
        Pointer to runtime.
    """

    # Load module from the toolkit folder.
    ds4sd_msg = load_tk_module(cmd_pointer, "DS4SD", "msgs", "ds4sd_msg")

    if GLOBAL_SETTINGS["display"] == "notebook":
        from tqdm.notebook import tqdm
    else:
        from tqdm import tqdm
    import pandas as pd

    cred_file = load_credentials(os.path.expanduser(f"{cmd_pointer.home_dir}/deepsearch_api.cred"))

    if cred_file["host"].strip() == "":
        host = DEFAULT_URL
    else:
        host = cred_file["host"]
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
    elif GLOBAL_SETTINGS["display"] != "notebook":
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
        output_error("No collection_key provided", return_val=False)
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
            "Domains": " / ".join(c.metadata.domain),
            "Collection Name": c.name,
            "Collection Key": c.source.index_key,
            "system_id": c.source.elastic_id,
        }
        for c in collections
    ]
    if val_elastic_id not in elastic_list:
        output_error("Invalid system_id, please choose from the following:", return_val=False)
        collectives = pd.DataFrame(result)
        output_table(collectives, is_data=False, return_val=False)
        return False
    if val_index_key not in index_list and val_index_key not in index_name_list:
        output_error("Invalid collection key or name, please choose from the following:", return_val=False)
        collectives = pd.DataFrame(result)
        output_table(collectives, is_data=False, return_val=False)
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
            coordinates=data_collection,  # The data collection to be queried
            aggregations=aggs,
        )
    else:
        query = DataQuery(
            search_query,  # The search query to be executed
            source=source_list,
            limit=page_size,  # The size of each request page
            coordinates=data_collection,  # The data collection to be queried
            aggregations=aggs,
        )

    # [Optional] Compute the number of total results matched.
    # This can be used to monitor the pagination progress.
    count_query = deepcopy(query)
    count_query.paginated_task.parameters["limit"] = 0
    count_results = api.queries.run(count_query)

    expected_total = count_results.outputs["data_count"]
    expected_pages = (expected_total + page_size - 1) // page_size  # This is simply a ceiling formula

    output_text("Estimated results: " + str(expected_total), return_val=False)
    if "estimate_only" in inputs:
        return None
    else:
        if expected_total > 100:
            if _confirm_prompt("Your query may take some time, do you wish to proceed?") is False:
                return None

    # Iterate through all results by fetching `page_size` results at the same time
    all_results = []
    all_aggs = {}
    try:
        cursor = api.queries.run_paginated_query(query)
        # raise Exception('This is a test error')
    except Exception as err:  # pylint: disable=broad-exception-caught
        output_error(ds4sd_msg("err_deepsearch", err), return_val=False)
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
        if all_aggs != {} and len(df.columns) > 1:
            output_text("<bold>Result distribution by year</bold>", pad_top=1, return_val=False)
            output_table(df.style.set_properties(**{"text-align": "left"}), is_data=False, return_val=False)
            print("")

    pd.set_option("display.max_colwidth", None)
    x = 0
    for row in all_results:
        while x < 20:
            x = x + 1
        result = {}
        if "_id" in row and GLOBAL_SETTINGS["display"] == "notebook" and "return_as_data" not in inputs:
            # result["ds_url"] = generate_url(host, data_collection, row["_id"])
            result["DS_URL"] = _make_clickable(_generate_url(host, data_collection, row["_id"]), "Deep Search Web Link")

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
                    result["arxivid"] = _make_clickable(f'https://arxiv.org/abs/{ref["value"]}', "ARXIVID Link")
                if ref["type"] == "doi":
                    result["doi"] = _make_clickable(f'https://doi.org/{ref["value"]}', "DOI Link")

        if edit_distance > 0:
            for field in row.get("highlight", {}).keys():
                for snippet in row["highlight"][field]:
                    result["Report"] = str(row["_source"]["file-info"]["filename"])
                    result["Field"] = field.split(".")[0]

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
        output_text("Search returned no result", return_val=False)
        return None
    if "save_as" in inputs:
        results_file = str(inputs["results_file"])
        if not results_file.endswith(".csv"):
            results_file = results_file + ".csv"

    df = pd.DataFrame(results_table)
    df = df.replace(np.nan, "", regex=True)
    if "save_as" in inputs:
        df.to_csv(
            cmd_pointer.workspace_path(cmd_pointer.settings["workspace"].upper()) + "/" + results_file, index=False
        )

    if GLOBAL_SETTINGS["display"] == "notebook":
        if "return_as_data" not in inputs:
            df = df.style
            df = df.set_properties(**{"text-align": "left"})
        else:
            return output_table(df, is_data=True).data

    elif GLOBAL_SETTINGS["display"] == "api":
        return df

    elif GLOBAL_SETTINGS["display"] == "terminal":
        if "save_as" not in inputs:
            df.style.format(hyperlinks="html")
            if "Title" in df:
                df["Title"] = df["Title"].str.wrap(50, break_long_words=True)
            if "Authors" in df:
                df["Authors"] = df["Authors"].str.wrap(25, break_long_words=True)
            if "Snippet" in df:
                df["Snippet"] = df["Snippet"].apply(lambda x: style(x))  # pylint diable:unnecessary-lambda=
                df["Snippet"] = df["Snippet"].str.wrap(70, break_long_words=True)
            if display_first > 0:
                df = df.truncate(after=display_first)

    return output_table(df, is_data=True)


def _confirm_prompt(question: str) -> bool:
    if GLOBAL_SETTINGS["display"] == "api":
        return True
    reply = None
    while reply not in ("y", "n"):
        reply = input(f"{question} (y/n): ").casefold()
        readline.remove_history_item(readline.get_current_history_length() - 1)
    return reply == "y"


def _make_clickable(url, name):
    if GLOBAL_SETTINGS["display"] == "notebook":
        return f'<a href="{url}"  target="_blank"> {name} </a>'
    else:
        return url


def _generate_url(host, data_source, document_hash, item_index=None):
    if isinstance(data_source, ElasticProjectDataCollectionSource):
        proj_key = data_source.proj_key
        index_key = data_source.index_key
        select_coords = {
            "privateCollection": index_key,
        }
        url = f"{host}/projects/{proj_key}/library/private/{index_key}"
    elif isinstance(data_source, ElasticDataCollectionSource):
        # TODO: remove hardcoding of community project
        proj_key = "1234567890abcdefghijklmnopqrstvwyz123456"
        index_key = data_source.index_key
        select_coords = {
            "collections": [index_key],
        }
        url = f"{host}/projects/{proj_key}/library/public"

    hash_expr = f'file-info.document-hash: "{document_hash}"'
    search_query = {
        **select_coords,
        "type": "Document",
        "expression": hash_expr,
        "filters": [],
        "select": [
            "_name",
            "description.collection",
            "prov",
            "description.title",
            "description.publication_date",
            "description.url_refs",
        ],
        "itemIndex": 0,
        "pageSize": 10,
        "searchAfterHistory": [],
        "viewType": "snippets",
        "recordSelection": {
            "record": {
                "id": document_hash,
            },
        },
    }
    if item_index is not None:
        search_query["recordSelection"]["itemIndex"] = item_index

    encoded_query = urllib.parse.quote(
        base64.b64encode(urllib.parse.quote(json.dumps(search_query, separators=(",", ":"))).encode("utf8")).decode(
            "utf8"
        )
    )

    url = f"{url}?search={encoded_query}"

    return url
