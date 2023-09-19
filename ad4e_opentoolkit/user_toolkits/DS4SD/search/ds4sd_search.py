import deepsearch as ds
from copy import deepcopy
from numerize import numerize
from deepsearch.cps.client.components.elastic import ElasticDataCollectionSource
from deepsearch.cps.queries import DataQuery
from deepsearch.cps.client.components.queries import RunQueryError
from ad4e_opentoolkit.helpers.style_parser import style
from ad4e_opentoolkit.helpers.style_parser import style, print_s, strip_tags, tags_to_markdown, parse_tags
from  ad4e_opentoolkit.helpers.output import output_table as output_table
from  ad4e_opentoolkit.helpers.output import output_text as output_text
import shutil
import textwrap,re

_tableformat = 'simple'

aggs = {
    "by_year": {
        "date_histogram": {
            "field": "description.publication_date",
                "calendar_interval": "year",
                "format": "yyyy",
                "min_doc_count": 0
        }
    }
}

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
    highlight = {"fields": {"*": {}}}
    highlight["fragment_size"] = 0
    if 'return_as_data' in inputs:
        highlight["pre_tags"]=[""]
        highlight["post_tags"]=[""]
    elif cmd_pointer.notebook_mode==False:
        highlight["pre_tags"]=['<green>']
        highlight["post_tags"]=['</green>']
    else:
        highlight["pre_tags"] = ["<span style='font-weight: bold; background-color: #FFFF00'>"]
        highlight["post_tags"] = ["</span>"]

    edit_distance=0
    search_query = inputs["search_text"]

    val_index_key = inputs["index_key"][val]

    if 'elastic_id' in inputs:
        val_elastic_id = inputs["elastic_id"][val]
    if 'page_size' in inputs:
        page_size = inputs["page_size"][val]
    if 'edit_distance' in inputs:
        edit_distance= int(inputs["edit_distance"][val])
    else:
        edit_distance= 10

    data_collection = ElasticDataCollectionSource(elastic_id=val_elastic_id, index_key=val_index_key)

    api = cmd_pointer.login_settings['toolkits_api'][cmd_pointer.login_settings['toolkits'].index('DS4SD') ]
    # Prepare the data query
    source_list = []
    if edit_distance>0:
        search_query= search_query+' ~'+str(edit_distance)

    is_docs=False

    if "show_data" in inputs:
        for i in inputs["show_data"]:
            if i.upper() == "DATA":
                source_list.extend(["subject", "attributes", "identifiers"])
            if i.upper() == "DOCS":
                source_list.extend(["description.title", "description.authors","file-info.filename", "identifiers"])
                is_docs=True
    else:
        source_list = ["subject", "attributes", "identifiers","file-info.filename"]
    
    if is_docs ==False:
        edit_distance=0

    if edit_distance >0:
        query = DataQuery(
            search_query,  # The search query to be executed
            # source=["subject", "attributes", "identifiers"], # Which fields of documents we want to fetch
            source=source_list,
            limit=page_size,  # The size of each request page
            highlight=highlight,
            coordinates=data_collection,  # The data collection to be queries
            aggregations=aggs
        )
    else:
        query = DataQuery(
        search_query,  # The search query to be executed
        # source=["subject", "attributes", "identifiers"], # Which fields of documents we want to fetch
        source=source_list,
        limit=page_size,  # The size of each request page
        coordinates=data_collection, # The data collection to be queries
        aggregations=aggs
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
    all_aggs={}
    cursor = api.queries.run_paginated_query(query)

    for result_page in tqdm(cursor, total=expected_pages):

        all_results.extend(result_page.outputs["data_outputs"])
        for year in result_page.outputs["data_aggs"]["by_year"]["buckets"]:
                
                if  year['key_as_string'] not in all_aggs:
                    all_aggs[year['key_as_string']]=0
                
                all_aggs[year['key_as_string']]= all_aggs[year['key_as_string']]+int(year['doc_count'])
                
#df.plot.bar(y="doc_count", x="key", figsize=(15, 5), xlabel="Language", ylabel="Number of reports", rot=0, legend=False, title="Number of non-English reports by language")

    results_table = []
    result = None
    import pandas as pd
    if is_docs:
        if cmd_pointer.notebook_mode==True:
            import IPython.display as display
            df= pd.json_normalize(all_aggs)
            if  len(df.columns)>1:
              
                display.display(output_text("<h3>Distribution of Returned Documents by Year</h3>"))
               
                display.display(output_table( df,cmd_pointer=cmd_pointer))
        elif all_aggs !={}:
            
            
            df= pd.json_normalize(all_aggs)
            if  len(df.columns)>1:
               
                output_text("\n<on_green>Distribution of Returned Documents by Year</on_green>")
        
                output_table( df,cmd_pointer=cmd_pointer)

    #df.plot.hist(ylabel="Number of reports",legend=False,figsize=(15, 5), title="Document Age Summary")
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
        if edit_distance>0:
            for field in row.get("highlight", {}).keys() :
                for snippet in row["highlight"][field]:
                    result["Snippet"]= re.sub(' +', ' ',snippet)
                    #result["Report"]= str(row["_source"]["file-info"]["filename"])
                    #result["Field"]= field
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
        if "identifiers" in row["_source"]:
            for ref in row["_source"]["identifiers"]:
                if ref["type"] == "arxivid":
                    result["arxivid"] = f'https://arxiv.org/abs/{ref["value"]}'
                if ref["type"] == "doi":
                    result["doi"] = f'https://doi.org/{ref["value"]}'
        if edit_distance>0:
            for field in row.get("highlight", {}).keys() :
                for snippet in row["highlight"][field]:
                    #result["Snippet"]= re.sub(' +', ' ',snippet)
                    result["Report"]= str(row["_source"]["file-info"]["filename"])
                    result["Field"]= field           
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
    if 'save_as' in inputs:
            results_file = str(inputs['results_file'])
            if not results_file.endswith('.csv'):
                results_file=results_file+'.csv'
    if cmd_pointer.notebook_mode == True:
        import pandas as pd
        df = pd.DataFrame(results_table)
      
        import numpy as np
        if 'save_as' in inputs:
            df.to_csv(cmd_pointer.workspace_path(cmd_pointer.settings['workspace'].upper()) + "/" + results_file,index=False)
        df =df.replace(np.nan, '', regex=True)
        
        
        if "return_as_data" in inputs:
            return (df)
        else:
            df =df.style.format(hyperlinks='html')
            df = df.set_properties(**{"text-align": "left"})
            return df

    else:
    
        import pandas as pd
        import numpy as np
        cmd_line_result = pd.DataFrame(results_table)

        if 'save_as' in inputs:

            cmd_line_result.to_csv(
                cmd_pointer.workspace_path(
                    cmd_pointer.settings['workspace'].upper()) +
                "/" + results_file,
                index=False)
        else:
            
            cmd_line_result.style.format(hyperlinks='html')
            if 'Title' in cmd_line_result:
                cmd_line_result['Title'] = cmd_line_result['Title'].str.wrap(50,break_long_words=True)
            if 'Authors' in cmd_line_result:
                cmd_line_result['Authors'] = cmd_line_result['Authors'].str.wrap(25,break_long_words=True)
            if 'Snippet' in cmd_line_result:
                cmd_line_result['Snippet']=cmd_line_result['Snippet'].apply(lambda x: style(x))
                
                cmd_line_result['Snippet'] = cmd_line_result['Snippet'].str.wrap(70,break_long_words=True)
                
                
                
                
        
        
            output_table(cmd_line_result.replace(np.nan, '', regex=True),cmd_pointer,tablefmt=_tableformat   )



def confirm_prompt(question: str) -> bool:
    import readline
    reply = None
    while reply not in ("y", "n"):
        reply = input(f"{question} (y/n): ").casefold()
        readline.remove_history_item(readline.get_current_history_length() - 1)
    return (reply == "y")
