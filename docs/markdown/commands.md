---
title: Commands
layout: home
nav_order: 4
---

<!--

DO NOT EDIT
-----------
This file auto-generated.
To update it, see openad/docs/generate-docs.py

-->

### Table of Contents
- [OpenAD](#openad)
  - [General](#general)
  - [Workspaces](#workspaces)
  - [Toolkits](#toolkits)
  - [Runs](#runs)
  - [Utility](#utility)
  - [LLM](#llm)
  - [File System](#file-system)
  - [Help](#help)
- [DS4SD](#ds4sd)
  - [Search Collection](#search-collection)
  - [Collections](#collections)
  - [Search Molecules](#search-molecules)
  - [Search Collections](#search-collections)
- [RXN](#rxn)
  - [Queries](#queries)
  - [Prediction](#prediction)
  - [Retrosynthesis](#retrosynthesis)
- [ST4SD](#st4sd)
- [GT4SD](#gt4sd)
  - [Generative Toolkit](#generative-toolkit)

## OpenAD

<details markdown="block">
<summary>See commands</summary>

### General

`openad`{: .cmd }
Display the openad splash screen.<br><br>

`get status`{: .cmd }
Display the currently selected workspace and toolkit.<br><br>

`display history`{: .cmd }
Display the last 30 commands run in your current workspace.<br><br>

`clear sessions`{: .cmd }
Clear any other sessions that may be running.<br><br>

<br>

### Workspaces

`set workspace <workspace_name>`{: .cmd }
Change the current workspace.<br><br>

`get workspace [ <workspace_name> ]`{: .cmd }
Display details a workspace. When no workspace name is passed, details of your current workspace are displayed.<br><br>

`create workspace <workspace_name> [ description('<description>') on path '<path>' ]`{: .cmd }
Create a new workspace with an optional description and path.<br><br>

`remove workspace <workspace_name>`{: .cmd }
Remove a workspace from your registry. Note that this doesn't remove the workspace's directory.<br><br>

`list workspaces`{: .cmd }
Lists all your workspaces.<br><br>

<br>

### Toolkits

`ds4sd`{: .cmd }
Display the splash screen for the DS4SD toolkit.<br><br>

`rxn`{: .cmd }
Display the splash screen for the RXN toolkit.<br><br>

`st4sd`{: .cmd }
Display the splash screen for the ST4SD toolkit.<br><br>

`gt4sd`{: .cmd }
Display the splash screen for the GT4SD toolkit.<br><br>

`list toolkits`{: .cmd }
List all installed toolkits. To see all available toolkits, run `list all toolkits`.<br><br>

`list all toolkits`{: .cmd }
List all available toolkits.<br><br>

`add toolkit <toolkit_name>`{: .cmd }
Install a toolkit.<br><br>

`remove toolkit <toolkit_name>`{: .cmd }
Remove a toolkit from the registry.
Note: This doesn't delete the toolkit code. If the toolkit is added again, a backup of the previous install is created in the toolkit directory at ~/.openad/toolkits.<br><br>

`set context <toolkit_name> [ reset ]`{: .cmd }
Set your context to the chosen toolkit. By setting the context, the selected toolkit functions become available to you. The optional parameter 'reset' can be used to reset your login information.<br><br>

`get context`{: .cmd }
Display the currently selected toolkit.<br><br>

`unset context`{: .cmd }
Exit your toolkit context. You will no longer have access to toolkit-specific functions.<br><br>

<br>

### Runs

`create run`{: .cmd }
Start recording a run.<br><br>

`save run as <run_name>`{: .cmd }
Stop recording a run and save it.<br><br>

`run <run_name>`{: .cmd }
Execute a previously recorded run. This will execute every command and continue regardless of any failures.<br><br>

`list runs`{: .cmd }
List all runs saved in the current workspace.<br><br>

`display run <run_name>`{: .cmd }
Display the commands stored in a certain run.<br><br>

<br>

### Utility

`display data '<csv_filename>'`{: .cmd }
Display data from a csv file.<br><br>

`-> result save [as '<csv_filename>']`{: .cmd }
Save table data to csv file.<br><br>

`-> result open`{: .cmd }
Explore table data in the browser.<br><br>

`-> result edit`{: .cmd }
Edit table data in the browser.<br><br>

`-> result copy`{: .cmd }
Copy table data to clipboard, formatted for spreadheet.<br><br>

`-> result display`{: .cmd }
Display the result in the CLI.<br><br>

`edit config '<json_config_file>' [ schema '<schema_file>']`{: .cmd }
Edit any JSON file in your workspace directly from the CLI. If a schema is specified, it will be used for validation and documentation.<br><br>

`show molecules using ( file '<mols_file>' | dataframe <dataframe> )
    [ save as '<sdf_or_csv_file>' | as molsobject ]`{: .cmd }
Launch the molecule viewer to examine and select molecules from a SMILES sdf/csv dataset.

Examples:

`show molecules using file 'base_molecules.sdf' as molsobject`
`show molecules using dataframe my_dataframe save as 'selection.sdf'`<br><br>

<br>

### LLM

`tell me <how to do xyz>`{: .cmd }
Ask your AI assistant how to do anything in OpenAD.<br><br>

`set llm  <language_model_name>`{: .cmd }
Set the target language model name for the "tell me" command.<br><br>

`clear llm auth`{: .cmd }
Clear the language model's authentication file.<br><br>

<br>

### File System

`list files`{: .cmd }
List all files in your current workspace.<br><br>

`import from '<external_source_file>' to '<workspace_file>'`{: .cmd }
Import a file from outside OpenAD into your current workspace.<br><br>

`export from '<workspace_file>' to '<external_file>'`{: .cmd }
Export a file from your current workspace to anywhere on your hard drive.<br><br>

`copy file '<workspace_file>' to '<other_workspace_name>'`{: .cmd }
Export a file from your current workspace to another workspace.<br><br>

`remove '<filename>'`{: .cmd }
Remove a file from your current workspace.<br><br>

<br>

### Help

`intro`{: .cmd }
Display an introduction to the OpenAD CLI.<br><br>

`docs`{: .cmd }
Open the documentation webpage.<br><br>

`?`{: .cmd }
List all available commands.<br><br>

`<soft>...</soft> ?`{: .cmd }
Display what a command does, or list all commands that contain this string.<br><br>

<br>

</details>

## DS4SD


<details markdown="block">
<summary>See commands</summary>

### Search Collection

`search collection '<collection name or key>' for '<search string>' using ( [ page_size=<int> system_id=<system_id> edit_distance=<integer> display_first=<integer>]) show (data|docs) [ estimate only|return as data|save as '<csv_filename>' ]`{: .cmd }
external

NOTE: The Using Clause Requires all the Parameters added to the Using Clause be in the defined order as per in the above help documentation<br><br>

`search collection '<collection name or key>' for '<search string>' using ( [ page_size=<int> system_id=<system_id> edit_distance=<integer> display_first=<integer>]) show (data|docs) [ estimate only|return as data|save as '<csv_filename>' ]`{: .cmd }
Performs a document search of the Deep Search repository based on a given collection. The required USING clause specifies the collection to search. Use 'estimate only' to perform a general search, returning the potential number of hits.

Parameters:
- `'<collection name or key>'` the name or index key for a collection. Use the command `display all collections` to identify collections.
- `'<search string>'` the search string for the search.

Search String Syntax: DeepSearch uses Elastic Search string query syntax, supporting operators like the following:
-- `+` signifies AND operation
-- `|` signifies OR operation
-- `-` negates a single token
-- `"` wraps a number of tokens to signify a phrase for searching
-- `*` at the end of a term signifies a prefix query
-- `(` and `)` signify precedence
-- `~N` after a word signifies edit distance (fuzziness)
-- `~N` after a phrase signifies slop amount

`USING` clause Options:
- `page_size=`<integer>`` - result pagination, the default is None.
- `system_id=`<system_id>` ` - system cluster id, the default is the value 'default'.
- `edit_distance=`<integer>``  - Set the search word span criteria for key words for document searches. the default is 5. When set to 0, no snippets will be be returned.
- `display_first=`<integer>`` - If display_first > 0, the displayed result set will be truncated at the given number. The default is 0.

Clauses:
- `show (data | docs ) ` - `data` Display structured data from within the documents or `docs` Display document context.
It is permitted to specify both in a single command e.g. ` show (data docs)`
- `estimate only` - Determine the potential number of hits.
- `return as data` - For Notebook or API mode. Removes all styling from the Pandas DataFrame, ready for further processing.

Examples:
- Look for documents that contain discussions on power conversion efficiency:

`search collection 'arxiv-abstract' for 'ide("power conversion efficiency" OR PCE) AND organ*' using ( edit_distance=20 system_id=default) show (docs)`

- Search the pubchem archive for 'Ibuprofen' and display related molecules' data:

`search collection 'pubchem' for 'Ibuprofen' show (data)`

- Search for patents which mention a specific smiles molecule:

`search collection 'uspto-patent' for 'identifiers._name:"smiles#ccc(coc(=o)cs)(c(=o)c(=o)cs)c(=o)c(=o)cs"' show (data)`


NOTE: The Using Clause Requires all the Parameters added to the Using Clause be in the defined order as per in the above help documentation<br><br>

<br>

### Collections

`display all collections [save as '<csv_file_name>']`{: .cmd }
This function displays all available collections in Deep Search.
If you use the `SAVE AS` clause, it will save a csv file to the current workspace.<br><br>

`display collections in domains from list [<list_of_domains>] [save as '<csv_file_name>']`{: .cmd }
This function displays collections that belong to the listed domains.
If you use the `SAVE AS` clause, it will save a csv file to the current workspace.<br><br>

`display collection details '<collection_name>' | '<collection_key>'`{: .cmd }
This function displays the details for a specified collection. You can specify either the name of a collection `<collection_name>` or its index key `<collection_key>`.<br><br>

`display collections for domain '<domain_name>'`{: .cmd }
This command displays the available collections in a given Deep Search `<domain_name>`.<br><br>

<br>

### Search Molecules

`search for similar molecules to '<smiles_string>' [save as '<csv_file_name>']`{: .cmd }
This command searches for molecules that are similar to the provided molecule or molecule substructure `<smiles_string>` provided.

For example `search for similar molecules to 'C1(C(=C)C([O-])C1C)=O'`

If you use the `SAVE AS` clause, it will save a csv file to the current workspace.<br><br>

`search for patents containing molecule ['<smiles_molecule>'| '<inchi_molecule>'] [save as '<csv_file_name>']`{: .cmd }
This command searches for mentions of a specified molecules in registered patents.
As input parameters you can provide either a SMILES version of a molecule `<smiles_molecule>` or Inchi `<inchi_molecule>`, which can either be in key or string format.

` search for patents containing molecule 'CC(C)(c1ccccn1)C(CC(=O)O)Nc1nc(-c2c[nH]c3ncc(Cl)cc23)c(C#N)cc1F' `

If you use the `SAVE AS` clause, it will save a csv file to the current workspace.<br><br>

`search for molecules in patents from [list ['<patent1>', '<patent2>' .....] | dataframe <dataframe_name> | file '<workspace_file name>'] [save as '<csv_file_name>']`{: .cmd }
This command searches for molecules that are mentioned in the defined list of patents. If sourcing patents are from CSV or dataframe, these must contain a column with 'PATENT ID' or 'patent id' as the heading.

For Example: ` search for molecules in patents from list ['CN108473493B','US20190023713A1'] `

If you use the `SAVE AS` clause, it will save a csv file to the current workspace.<br><br>

`search for substructure instances of '<smiles_string>' [save as '<csv_file_name>']`{: .cmd }
This command searches for molecules with the instance of a molecule in their substructure, as defined in the `<smiles_string>` string.
If you use the `SAVE AS` clause, it will save a csv file to the current workspace.

For example: ` search for substructure instances of 'C1(C(=C)C([O-])C1C)=O' save as 'my_mol'`<br><br>

<br>

### Search Collections

`display collection matches for '<search_string>' [save as '<csv_file_name>']`{: .cmd }
This command searches all collections for documents that contain a given Deep Search `<search_string>`. This helps choose document collection(s) for subsequent search. Use `<index_key>` from the returned table in a search.
If you use the `SAVE AS` clause, it will save a csv file to the current workspace.<br><br>

<br>

</details>

## RXN


<details markdown="block">
<summary>See commands</summary>

### Queries

`list rxn models`{: .cmd }
lists current rxn AI Models available to the user<br><br>

<br>

### Prediction

`predict reaction topn in batch from (dataframe <dataframe_name> | file '<file_name.csv>' | list ['#smiles_reaction','#smiles_reaction') [USING (topn=<integer> ai_model='<existing_model>')] [use_saved]`{: .cmd }
This command performs a reaction prediction for topn providing results for a given list of reactions. The list of reactions can be specified as a string list, data frame or csv file from the current workspace. For data frames and csv files it will take the column with the name ‘reactions’.

In the `FROM` clause reactions are defined by a list of reactions where are SMILES string is delimited by '.' e.g. `'BrBr.c1ccc2cc3ccccc3cc2c1'`

The optional `USING` clause can specify an AI model, a value for topn, or both:
- `ai_model=’`<model_name>`’ ` The default value is '2020-07-01'
- `topn=`<integer>``  this sets the top n results, the default value is 3

Examples:
`predict reaction topn batch from list ['BrBr.c1ccc2cc3ccccc3cc2c1CCO' , 'BrBr.c1ccc2cc3ccccc3cc2c1']`

`predict reaction topn batch from list ['BrBr.c1ccc2cc3ccccc3cc2c1CCO' , 'BrBr.c1ccc2cc3ccccc3cc2c1'] using ( topn=6)`

You can also use previously generated results buy optionally using `use_saved` at the end of the command and it will use the results of any previously run commands with the same parameters while the toolkit has been installed.

`predict reaction topn batch from list ['BrBr.c1ccc2cc3ccccc3cc2c1CCO' , 'BrBr.c1ccc2cc3ccccc3cc2c1'] using (topn=6) use_saved `


NOTE: The Using Clause Requires all the Parameters added to the Using Clause be in the defined order as per in the above help documentation<br><br>

`predict reaction '<reaction-smiles-string>' [USING (ai_model='<valid_ai_model>')] [use_saved]`{: .cmd }
This command 'forward predicts' a reaction for a given SMILES string.

In the `FROM` clause is a list of reactions: SMILES strings delimited by a period '.', e.g. `'BrBr.c1ccc2cc3ccccc3cc2c1'`


The optional `USING` clause specifies a particular AI model.
-`ai_model=’`<model_name>`’` The default value is '2020-07-01'

Example:
`predict reaction 'BrBr.c1ccc2cc3ccccc3cc2c1CCO'`

You can optionally use previously generated results with `use_saved` at the end of the command. It will use the results of any previous commands run with the same parameters.

`predict reaction 'BrBr.c1ccc2cc3ccccc3cc2c1CCO' use_saved`

NOTE: The Using Clause Requires all the Parameters added to the Using Clause be in the defined order as per in the above help documentation<br><br>

`predict reaction in batch from [dataframe < dataframe_name > ] | [file '<file_name.csv>'] | [list ['#smiles','#smiles']]  [USING ( ai_model='<existing_model>')] [use_saved]`{: .cmd }
This command performs a reaction prediction providing results for a given list of possible reaction paths. The list of reactions can be specified as a string list, data frame or csv file from the current workspace. For data frames and csv files it will take the column with the name 'reactions'.

In the `FROM` clause reactions are defined by a list of reactions where are SMILES string is delimited by '.' e.g. `'BrBr.c1ccc2cc3ccccc3cc2c1'`

The optional `USING` clause specifies an AI model other than the default model.
- `ai_model=’`<model_name>`’ `The default ai_model is '2020-07-01'
Examples:
`predict reaction batch from list ['BrBr.c1ccc2cc3ccccc3cc2c1CCO' , 'BrBr.c1ccc2cc3ccccc3cc2c1']`

You can also use previously generated results by optionally using `use_saved` at the end of the command and it will use the results of any previously run commands with the same parameters while the toolkit has been installed.

`predict reaction batch from list ['BrBr.c1ccc2cc3ccccc3cc2c1CCO' , 'BrBr.c1ccc2cc3ccccc3cc2c1'] use_saved`

NOTE: The Using Clause Requires all the Parameters added to the Using Clause be in the defined order as per in the above help documentation<br><br>

<br>

### Retrosynthesis

`interpret recipe '<recipe-paragraph> | <workspace-file>'`{: .cmd }
This command builds a set of actions interpreted from a provided recipe defined as a provided string or a file in the current workspace in the parameter ``<recipe-paragraph>` | `<workspace-file>``<br><br>

`predict retrosynthesis '<product_SMILES_string>' [USING ( option=<valid_input> option2=<valid_input> )]`{: .cmd }
This command performs automatic retro synthesis route prediction on a given molecule.

The parameter `'<product_SMILES_string>'` takes a valid SMILES string.

Options for `USING` clause are:
- `availability_pricing_threshold=`<int>` ` maximum price in USD per g/ml of compounds. Default: no threshold.
- `available_smiles='<list of SMILES>'` list of molecules available as precursors, with delimiter '.'
- `exclude_smiles='<list of SMILES>'` list of molecules to exclude from the set of precursors, delimiter '.'
- `exclude_substructures='<list of SMILES>'` substructures to exclude, delimiter '.'
- `exclude_target_molecule=`<boolean>`` excluded target molecule, default True
- `fap=`<float>`` Every retrosynthetic step is evaluated with the FAP, a step is retained when forward confidence is greater than FAP, default 0.6
- `max_steps=`<int>`` The max steps, default is 3
- `nbeams=`<int>` ` The maximum number of beams exploring the hypertree, default 10
- `pruning_steps=`<int>`` The number of steps to prune a hypertree, default 2
- `ai_model='<ai_model_name>'` default '2020-07-01'

An example command is:
`predict retrosynthesis 'BrCCc1cccc2c(Br)c3ccccc3cc12' using (max_steps=3) `

NOTE: The Using Clause Requires all the Parameters added to the Using Clause be in the defined order as per in the above help documentation<br><br>

<br>

</details>

## ST4SD


<details markdown="block">
<summary>See commands</summary>

</details>

## GT4SD


<details markdown="block">
<summary>See commands</summary>

### Generative Toolkit

`exec inference()`{: .cmd }
this is a gt4sd test function<br><br>

<br>

</details>
