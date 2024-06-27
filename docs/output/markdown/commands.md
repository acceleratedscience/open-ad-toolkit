---
title: Commands
layout: home
nav_order: 4
---

<!--

DO NOT EDIT
-----------
This file auto-generated.
To update it, see openad/docs/generate_docs.py

-->

### Table of Contents
- [OpenAD](#openad)
  - [General](#general)
  - [Workspaces](#workspaces)
  - [Molecules](#molecules)
  - [Utility](#utility)
  - [Toolkits](#toolkits)
  - [Runs](#runs)
  - [LLM](#llm)
  - [File System](#file-system)
  - [Help](#help)
  - [Model](#model)
- [DS4SD](#ds4sd)
  - [Search Molecules](#search-molecules)
  - [Search Collections](#search-collections)
  - [Collections](#collections)
- [RXN](#rxn)
  - [General](#general)
  - [Retrosynthesis](#retrosynthesis)
  - [Prediction](#prediction)
- [ST4SD](#st4sd)

## OpenAD

<details markdown="block">
<summary>See commands</summary>

### General

`openad`{: .cmd }
Display the openad splash screen. <br><br>

`get status`{: .cmd }
Display the currently selected workspace and toolkit. <br><br>

`display history`{: .cmd }
Display the last 30 commands run in your current workspace. <br><br>

`clear sessions`{: .cmd }
Clear any other sessions that may be running. <br><br>

<br>

### Workspaces

`set workspace <workspace_name>`{: .cmd }
Change the current workspace. <br><br>

`get workspace [ <workspace_name> ]`{: .cmd }
Display details a workspace. When no workspace name is passed, details of your current workspace are displayed. <br><br>

`create workspace <workspace_name> [ description('<description>') on path '<path>' ]`{: .cmd }
Create a new workspace with an optional description and path. <br><br>

`remove workspace <workspace_name>`{: .cmd }
Remove a workspace from your registry. Note that this doesn't remove the workspace's directory. <br><br>

`list workspaces`{: .cmd }
Lists all your workspaces. <br><br>

<br>

### Molecules

`add molecule <name> | <smiles> | <inchi> | <inchikey> | <cid> [ as '<name>' ] [ basic ] [force ]`{: .cmd }
This command is how you add a molecule to a current working list of molecules in memory. When adding a molecule by name, this name will become the molecule's identifying string.  <br> 

It will take any molecules identifier from the following categories: <br> 
    -`smiles ` <br> 
    -`name or synonym` <br> 
    -`smiles` <br> 
    -`inchi` <br> 
    -`inchikey ` <br> 
    -`cid ` <br> 

Options : <br> 
    - `as <name> `: if the ` as '<name>' ` not used the molecule the  molecule identfier will be used for the molecules name. if the ` as '<name>' ` not used the molecule the  molecule identfier will be used for the molecules name. <br> 
        You can set or override an name later for  any molecule with the `rename molecule` command. <br> 
    - ` basic ` Creates a molecule that does not have its properties and synonyms populated with pubchem data, this feature is only valid with a SMILES molecule identifier <br> 
    - `force`: The `force` option allows you to ovveride the confirmation that you wish to add a molecule. <br> 




You can use the 'mol' shorthand instead of 'molecule'. <br> 

You can specify any molecule by SMILES or InChI, and PubChem classified molecules also by name, InChIKey or their PubChem CID.  <br> 
 A molecule identifier can be in single quotes or defined with unquoted text. If you have spaces in your molecule identifier e.g. a name, then you must user a single quoted string <br> 

If you use the name of a molecule, the tool will do a caseless search of the names and synonyms first in current working list, then on pubchem. <br> 


Examples of how to add a molecule to your working list: <br> 
- Add a molecule by name: <br> 
`add molecule aspirin` <br> 
or with single quotes <br> 
` display molecule 'Aspirin 325 mg' ` <br> 

- Add a molecule by name and force through the acknowledgement to add it: <br> 
`add molecule aspirin force` <br> 

- Add a molecule by SMILES: <br> 
`add molecule CC(=O)OC1=CC=CC=C1C(=O)O` <br> 

- Add a molecule by SMILES without populated pubchem properties: <br> 
`add molecule CC(=O)OC1=CC=CC=C1C(=O)O basic ` <br> 

- Add a molecule by CID: <br> 
`add mol 2244` <br> 

- Add a molecule by InChI: <br> 
`add mol InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)` <br> 
 or with single quotes <br> 
 `add mol 'InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)'` <br> 

- Add a molecule by InChIKey: <br> 
`add mol BSYNRYMUTXBXSQ-UHFFFAOYSA-N` <br> 

- Add a molecule by InChIKey nd set its name to "mymol": <br> 
`add mol BSYNRYMUTXBXSQ-UHFFFAOYSA-N as 'mymol' ` <br> 

- Add a molecule by SMILES nd set its name to "mymol" and not prepopulate values from pubchem: <br> 
`add mol CC(=O)OC1=CC=CC=C1C(=O)O as 'mymol' basic ` <br><br>

`display molecule <name> | <smiles> | <inchi> | <inchikey> | <cid>`{: .cmd }
This command will display a molecule's identifiers, propoerties, synonyms and any Analysis results it has been enriched with. <br> 
if the molecule exists in the current molecule workling list in memory the molecule will be retrieved from there if not pubchem will be checked to see if the molecule and its information is avialable there. <br> 

You can use the 'mol' shorthand instead of 'molecule'. <br> 

If the requested molecule exists in your current working list, that version will be used. <br> 

You can specify any molecule by SMILES or InChI, and PubChem classified molecules also by name, InChIKey or their PubChem CID.  <br> 
 A molecule identifier can be in single quotes or defined with unquoted text. If you have spaces in your molecule identifier e.g. a name, then you must user a single quoted string <br> 

If you use the name of a molecule, the tool will do a caseless search of the names and synonyms first in current working list, then on pubchem. <br> 

Examples: <br> 
- Display a molecule by name: <br> 
`display molecule Aspirin` <br> 

- Display a molecule by SMILES: <br> 
`display molecule CC(=O)OC1=CC=CC=C1C(=O)O` <br> 

- Display a molecule by InChI: <br> 
`display mol InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)` <br> 

- Display a molecule by InChIKey string: <br> 
`display mol BSYNRYMUTXBXSQ-UHFFFAOYSA-N` <br> 

- Display a molecule by CID: <br> 
`display mol 2244` <br><br>

`display sources <name> | <smiles> | <inchi> | <inchikey> | <cid>`{: .cmd }
Display the sources of a molecule's properties, attributing back to how they were calculated or sourced. <br> 

If the requested molecule exists in your current working list, that version will be used. <br> 

You can specify any molecule by SMILES or InChI, and PubChem classified molecules also by name, InChIKey or their PubChem CID.  <br> 
 A molecule identifier can be in single quotes or defined with unquoted text. If you have spaces in your molecule identifier e.g. a name, then you must user a single quoted string <br> 

If you use the name of a molecule, the tool will do a caseless search of the names and synonyms first in current working list, then on pubchem. <br><br>

`rename molecule <molecule_identifer_string> as <molecule_name>`{: .cmd }
Rename a molecule in the current working list. <br> 

{MOL_SHORTHAND} <br> 

Example: <br> 
Let's say you've added a molecule "CC(=O)OC1=CC=CC=C1C(=O)O" to your current working list of molecules, you can then rename it as such: <br> 
`rename molecule CC(=O)OC1=CC=CC=C1C(=O)O as Aspirin` <br><br>

`export molecule <name> | <smiles> | <inchi> | <inchikey> | <cid> [ as file ]`{: .cmd }
When run inside a jupyter lab notebook, this will return a dictionary of the molecule's properties. When run from the command line, or when `as file` is set, the molecule will be saved to your workspace as a JSON file, named after the molecule's identifying string. <br> 
If the molecule is in your current working list it will be retrieved from there, if the molecule is not there pubchem will be called to retrieve the molecule. <br> 

You can use the 'mol' shorthand instead of 'molecule'. <br> 

If the requested molecule exists in your current working list, that version will be used. <br> 

If you use the name of a molecule, the tool will do a caseless search of the names and synonyms first in current working list, then on pubchem. <br> 

Examples <br> 
- `export molecule aspirin` <br> 
- `export molecule aspirin as file` <br><br>

`remove molecule <name> | <smiles> | <inchi> | <inchikey> | <cid>`{: .cmd }
Remove a molecule from the current working list based on a given molecule identifier. <br> 

{MOL_SHORTHAND} <br> 

Examples: <br> 
- Remove a molecule by name: <br> 
`remove molecule Aspirin` <br> 

- Remove a molecule by SMILES: <br> 
`remove molecule CC(=O)OC1=CC=CC=C1C(=O)O` <br> 

- Remove a molecule by InChIKey: <br> 
`remove mol  BSYNRYMUTXBXSQ-UHFFFAOYSA-N` <br> 

- Remove a molecule by InChI: <br> 
`remove mol  InChI='1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)'` <br> 

- Remove a molecule by CID: <br> 
`remove mol 2244` <br><br>

`list molecules`{: .cmd }
List all molecules in the current working list. <br><br>

`save molecule-set as <molecule_set_name>`{: .cmd }
Save the current molecule workking list to a molecule-set in your workspace. <br> 

Example: <br> 
`save molset as my_working_set` <br><br>

`load molecule-set|molset <molecule-set_name>`{: .cmd }
Loads a molecule-set from your workspace, and replaces your current list of molecules with the molecules from the given  molecule-set. <br> 
Example: <br> 
`load molset my_working_set` <br><br>

`merge molecule-set|molset <molecule-set_name> [merge only] [append only]`{: .cmd }
This command merges a molecule-set from your workspace into cour current working list of molecules in memory, and updates properties/Analysis in existing molecules or appends new molecules to the working list. <br> 

Options: <br> 
    - ` merge only` Only merges with existing molecules in list <br> 
    - ` append only` Only append molecules not in list <br> 
`merge molset my_working_set` <br><br>

`list molecule-sets`{: .cmd }
List all molecule sets in your workspace. <br><br>

`enrich molecules with analysis`{: .cmd }
This command Enriches every molecule in your current working list of molecules with the analysis results. This assumes that molecules in the current working list was the input or result for the analysis. <br> 

            This command currently covers results from the following Analysis commands: <br> 
            - RXN Toolkit `predict Reaction` <br> 
            - RXN Toolkit `predict retrosynthesis ` <br> 
            - DS4SD Toolkit `search for patents containing molecule` <br> 
            - DS4SD Toolkit `search for similiar molecules` <br> 

            See the Deep Search toolkit  and RXN toolkit help for further assistance on these commands.  <br><br>

`clear analysis cache`{: .cmd }
this command clears the cache of analysis results for your current workspace. <br><br>

`clear molecules`{: .cmd }
This command clears the working list of molecules. <br><br>

`@(<name> | <smiles> | <inchi> | <inchikey> | <cid>)>><molecule_property_name>`{: .cmd }
This command request the given property of a molecule, it will first try and retrieve the provided molecule from your working list of molecules, if it is not there it will will try and retrieve the molecule from pubchem. <br> 

The `@` symbol should be followed by the molecule's name, SMILES, InChI, InChIKey or CID, then after the `>>` include one of the properties mentioned below. <br> 

E.g. `@aspirin>>xlogp` <br> 

You can specify any molecule by SMILES or InChI, and PubChem classified molecules also by name, InChIKey or their PubChem CID.  <br> 
 A molecule identifier can be in single quotes or defined with unquoted text. If you have spaces in your molecule identifier e.g. a name, then you must user a single quoted string <br> 

If you use the name of a molecule, the tool will do a caseless search of the names and synonyms first in current working list, then on pubchem. <br> 

Examples of how to retrieve the value of a molecules property: <br> 
- Obtain the molecular weight of the molecule known as Aspirin. <br> 
`@aspirin>>molecular_weight` <br> 

- Obtain the canonical smiles string for a molecule known as Aspirin. <br> 
`@aspirin>>canonical_smiles` <br> 

- Obtain a molecules xlogp value using a SMILES string. <br> 
`@CC(=O)OC1=CC=CC=C1C(=O)O>>xlogp` <br> 

Available properties: `atom_stereo_count`, `bond_stereo_count`, `canonical_smiles`, `charge`, `cid`, `complexity`, `conformer_count_3d`, `conformer_id_3d`, `conformer_model_rmsd_3d`, `conformer_rmsd_3d`, `coordinate_type`, `covalent_unit_count`, `defined_atom_stereo_count`, `defined_bond_stereo_count`, `effective_rotor_count_3d`, `exact_mass`, `feature_acceptor_count_3d`, `feature_anion_count_3d`, `feature_cation_count_3d`, `feature_count_3d`, `feature_donor_count_3d`, `feature_hydrophobe_count_3d`, `feature_ring_count_3d`, `h_bond_acceptor_count`, `h_bond_donor_count`, `heavy_atom_count`, `inchi`, `inchikey`, `isomeric_smiles`, `isotope_atom_count`, `iupac_name`, `mmff94_energy_3d`, `mmff94_partial_charges_3d`, `molecular_formula`, `molecular_weight`, `molecular_weight_exact`, `monoisotopic_mass`, `multipoles_3d`, `multipoles_3d`, `pharmacophore_features_3d`, `pharmacophore_features_3d`, `rotatable_bond_count`, `sol`, `sol_classification`, `tpsa`, `undefined_atom_stereo_count`, `undefined_bond_stereo_count`, `volume_3d`, `x_steric_quadrupole_3d`, `xlogp`, `y_steric_quadrupole_3d`, `z_steric_quadrupole_3d` <br><br>

`load molecules using file '<csv_or_sdf_filename>' [ merge with pubchem ]`{: .cmd }
This command Loads molecules from a CSV or SDF file into the molecule working list. Optionally you can add `merge with pubchem` to the command to fill in missing properties of the molecule. <br><br>

`export molecules [ as <csv_filename> ]`{: .cmd }
This command exports the molecules in the current working list of molecules. <br> 

When run inside a Notebook, this will return a dataframe. When run from the command line, the molecules will be saved to your workspace as a CSV file named "result_#.csv". The rows will be numbered with the highest number representing the latest molecule that was added. <br><br>

`show molecule <name> | <smiles> | <inchi> | <inchikey> | <cid>`{: .cmd }
Inspect a molecule in the browser. <br> 

{MOL_SHORTHAND} <br> 

When you show a molecule by SMILES or InChI, we can display it immediately. When you show a molecule by name, InChIKey or PubChem CID, we need to first retrieve it from PubChem, which can take a few seconds. <br> 

Examples of how to show a molecule and its proerties in the molecule viewer: <br> 
- `show mol aspirin` <br> 
- `show mol CC(=O)OC1=CC=CC=C1C(=O)O` <br> 
- `show mol InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)` <br> 
- `show mol 2244` <br><br>

`show molecules using ( file '<mols_file>' | dataframe <dataframe> ) [ save as '<sdf_or_csv_file>' | as molsobject ]`{: .cmd }
Launch the molecule viewer to examine and select molecules from a SMILES sdf/csv dataset. <br> 

    Examples of how to show molecules in mols2grid: <br> 
    - `show molecules using file 'base_molecules.sdf' as molsobject` <br> 
    - `show molecules using dataframe my_dataframe save as 'selection.sdf'` <br><br>

<br>

### Utility

`load molecules using dataframe <dataframe> [ merge with pubchem ]`{: .cmd }
"             <br> 
This command Load molecules into the molecule working list from a dataframe.  <br> 

If the ` merge with pubchem`  clause is used then loaded molecules will have properties that are not in the source file filled in using pubchem requests, this will slow the process down <br><br>

`merge molecules data using dataframe <dataframe> [ merge with pubchem ]`{: .cmd }
"             <br> 
This command merges molecules into the molecule working list from a dataframe.  <br><br>

`display data '<filename.csv>'`{: .cmd }
Display data from a csv file. <br><br>

`-> result save [as '<filename.csv>']`{: .cmd }
Save table data to csv file. <br><br>

`-> result open`{: .cmd }
Explore table data in the browser. <br><br>

`-> result edit`{: .cmd }
Edit table data in the browser. <br><br>

`-> result copy`{: .cmd }
Copy table data to clipboard, formatted for spreadheet. <br><br>

`-> result display`{: .cmd }
Display the result in the CLI. <br><br>

`-> result as dataframe`{: .cmd }
Return the result as dataframe (only for Jupyter Notebook) <br><br>

`edit config '<json_config_file>' [ schema '<schema_file>']`{: .cmd }
Edit any JSON file in your workspace directly from the CLI. If a schema is specified, it will be used for validation and documentation. <br><br>

<br>

### Toolkits

`ds4sd`{: .cmd }
Display the splash screen for the DS4SD toolkit. <br><br>

`rxn`{: .cmd }
Display the splash screen for the RXN toolkit. <br><br>

`st4sd`{: .cmd }
Display the splash screen for the ST4SD toolkit. <br><br>

`list toolkits`{: .cmd }
List all installed toolkits. To see all available toolkits, run `list all toolkits`. <br><br>

`list all toolkits`{: .cmd }
List all available toolkits. <br><br>

`add toolkit <toolkit_name>`{: .cmd }
Install a toolkit. <br><br>

`remove toolkit <toolkit_name>`{: .cmd }
Remove a toolkit from the registry. <br> 

<b>Note:</b> This doesn't delete the toolkit code. If the toolkit is added again, a backup of the previous install is created in the toolkit directory at `~/.openad/toolkits`. <br><br>

`update toolkit <toolkit_name>`{: .cmd }
Update a toolkit with the latest version. It is recommended to do this on a regular basis. <br><br>

`update all toolkits`{: .cmd }
Update all installed toolkits with the latest version. Happens automatically whenever OpenAD is updated to a new version. <br><br>

`set context <toolkit_name> [ reset ]`{: .cmd }
Set your context to the chosen toolkit. By setting the context, the selected toolkit functions become available to you. The optional parameter `reset` can be used to reset your login information. <br><br>

`get context`{: .cmd }
Display the currently selected toolkit. <br><br>

`unset context`{: .cmd }
Exit your toolkit context. You will no longer have access to toolkit-specific functions. <br><br>

<br>

### Runs

`create run`{: .cmd }
Start recording a run. <br><br>

`save run as <run_name>`{: .cmd }
Stop recording a run and save it. <br><br>

`run <run_name>`{: .cmd }
Execute a previously recorded run. This will execute every command and continue regardless of any failures. <br><br>

`list runs`{: .cmd }
List all runs saved in the current workspace. <br><br>

`display run <run_name>`{: .cmd }
Display the commands stored in a certain run. <br><br>

<br>

### LLM

`tell me <how to do xyz>`{: .cmd }
Ask your AI assistant how to do anything in OpenAD. <br><br>

`set llm  <language_model_name>`{: .cmd }
Set the target language model name for the `tell me` command. <br><br>

`clear llm auth`{: .cmd }
Clear the language model's authentication file. <br><br>

<br>

### File System

`list files`{: .cmd }
List all files in your current workspace. <br><br>

`import from '<external_source_file>' to '<workspace_file>'`{: .cmd }
Import a file from outside OpenAD into your current workspace. <br><br>

`export from '<workspace_file>' to '<external_file>'`{: .cmd }
Export a file from your current workspace to anywhere on your hard drive. <br><br>

`copy file '<workspace_file>' to '<other_workspace_name>'`{: .cmd }
Export a file from your current workspace to another workspace. <br><br>

`remove '<filename>'`{: .cmd }
Remove a file from your current workspace. <br><br>

<br>

### Help

`intro`{: .cmd }
Display an introduction to the OpenAD CLI. <br><br>

`docs`{: .cmd }
Open the documentation webpage. <br><br>

`?`{: .cmd }
List all available commands. <br><br>

`? ...<soft>   --> List all commands containing "..."</soft>`{: .cmd }
<br>

`... ?<soft>   --> List all commands starting with "..."</soft>`{: .cmd }
<br>

<br>

### Model

`model auth list`{: .cmd }
show authentication group mapping <br><br>

`model auth add group '<auth_group>' with '<api_key>'`{: .cmd }
add an authentication group for model services to use <br><br>

`model auth remove group '<auth_group>'`{: .cmd }
remove an authentication group <br><br>

`model auth add service '<service_name>' to group '<auth_group>'`{: .cmd }
attach an authentication group to a model service <br><br>

`model auth remove service '<service_name>'`{: .cmd }
detatch an authentication group from a model service <br><br>

`model service status`{: .cmd }
get the status of currently cataloged services <br><br>

`model service describe '<service_name>'|<service_name>`{: .cmd }
get the configuration of a service <br><br>

`model catalog list`{: .cmd }
get the list of currently cataloged services <br><br>

`uncatalog model service '<service_name>'|<service_name>`{: .cmd }
uncatalog a model service  <br> 

 Example:  <br> 
`uncatalog model service 'gen'` <br><br>

`catalog model service from (remote) '<path or github>' as  '<service_name>'|<service_name>   USING (<parameter>=<value> <parameter>=<value>)`{: .cmd }
catalog a model service from a path or github or remotely from an existing OpenAD service. <br> 
(USING) optional headers parameters for communication with service backend. <br> 

Example: <br> 

-`catalog model service from 'git@github.com:acceleratedscience/generation_inference_service.git' as 'gen'` <br> 

or to catalog a remote service shared with you: <br> 
-`catalog model service from remote 'http://54.235.3.243:30001' as gen` <br><br>

`model service up '<service_name>'|<service_name> [no_gpu]}`{: .cmd }
launches a cataloged model service. <br> 
If you do not want to launch a service with GPU you should specify `no_gpu` at the end of the command. <br> 
Examples: <br> 

-`model service up gen` <br> 

-`model service up 'gen'` <br> 

-`model service up gen no_gpu` <br><br>

`model service local up '<service_name>'|<service_name>`{: .cmd }
launch a model service locally. <br> 

            Example: <br> 
              ` model service local up gen` <br><br>

`model service down '<service_name>'|<service_name>`{: .cmd }
Bring down a model service   <br> 
 Examples:  <br> 

`model service down gen`  <br> 

`model service down 'gen'`  <br><br>

<br>

</details>

## DS4SD


<details markdown="block">
<summary>See commands</summary>

### Search Molecules

`search for similar molecules to '<smiles>' [ save as '<filename.csv>' ]`{: .cmd }
Search for molecules that are similar to the provided molecule or molecule substructure as provided in the `<smiles_string>`. <br> 

Use the `save as` clause to save the results as a csv file in your current workspace. <br> 

Example: <br> 
`search for similar molecules to 'C1(C(=C)C([O-])C1C)=O'` <br><br>

`search for molecules in patents from list ['<patent1>', '<patent2>', ...] | dataframe <dataframe_name> | file '<filename.csv>' [ save as '<filename.csv>' ]`{: .cmd }
Search for molecules mentioned in a defined list of patents. When sourcing patents from a CSV or DataFrame, there must be column named "PATENT ID" or "patent id". <br> 

Use the `save as` clause to save the results as a csv file in your current workspace. <br> 

Example: <br> 
`search for molecules in patents from list ['CN108473493B','US20190023713A1']` <br><br>

`search for patents containing molecule '<smiles>' | '<inchi>' | '<inchikey>' [ save as '<filename.csv>' ]`{: .cmd }
Search for mentions of a specified molecules in registered patents. The queried molecule can be described as a SMILES string, InChI or InChiKey. <br> 

Use the `save as` clause to save the results as a csv file in your current workspace. <br> 

Example: <br> 
`search for patents containing molecule 'CC(C)(c1ccccn1)C(CC(=O)O)Nc1nc(-c2c[nH]c3ncc(Cl)cc23)c(C#N)cc1F'` <br><br>

`search for substructure instances of '<smiles>' [ save as '<filename.csv>' ]`{: .cmd }
Search for molecules by substructure, as defined by the `<smiles_string>`. <br> 

Use the `save as` clause to save the results as a csv file in your current workspace. <br> 

Example: <br> 
`search for substructure instances of 'C1(C(=C)C([O-])C1C)=O' save as 'my_mol'` <br><br>

<br>

### Search Collections

`search collection '<collection_name_or_key>' for '<search_string>' [ using (page_size=<int> system_id=<system_id> edit_distance=<integer> display_first=<integer>) ] show (data | docs) [ estimate only | return as data | save as '<filename.csv>' ]`{: .cmd }
Performs a document search of the Deep Search repository based on a given collection. The required `using` clause specifies the collection to search. Use `estimate only` to return only the potential number of hits. <br> 

Parameters: <br> 
- `<collection_name_or_key>` The name or index key for a collection. Use the command `display all collections` to list available collections. <br> 
- `<search_string>` The search string for the search. <br> 

The `<search_string>` supports elastic search string query syntax: <br> 
- `+` Signifies AND operation. <br> 
- `|` Signifies OR operation. <br> 
- `-` Negates a single token. <br> 
- `\"` Wraps a number of tokens to signify a phrase for searching. <br> 
- `*` At the end of a term -> signifies a prefix query <br> 
- `(` & `)` Signifies precedence <br> 
- `~N` After a word -> signifies edit distance (fuzziness) <br> 
- `~N` After a phrase -> signifies slop amount <br> 

Options for the `using` clause: <br> 
  > **Note:** The `using` clause requires all enclosed parameters to be defined in the same order as listed below. <br> 

- `page_size=<integer>` Result pagination, the default is None. <br> 
- `system_id=<system_id>` System cluster id, the default is 'default'. <br> 
- `edit_distance=<integer>` (0-5) Sets the search word span criteria for key words for document searches, the default is 5. When set to 0, no snippets will be be returned. <br> 
- `display_first=<integer>` When set, the displayed result set will be truncated at the given number. <br> 

Clauses: <br> 
- `show (data | docs)`: <br> 
    - `data` Display structured data from within the documents. <br> 
    - `docs` Display document context and preview snippet. <br> 
    Both can be combined in a single command, e.g. `show (data docs)` <br> 
- `estimate only` Determine the potential number of hits. <br> 
- `return as data` For Notebook or API mode. Removes all styling from the Pandas DataFrame, ready for further processing. <br> 

Examples: <br> 
- Look for documents that contain discussions on power conversion efficiency: <br> 
`search collection 'arxiv-abstract' for 'ide(\"power conversion efficiency\" OR PCE) AND organ*' using ( edit_distance=20 system_id=default) show (docs)` <br> 

- Search the PubChem archive for 'Ibuprofen' and display related molecules' data: <br> 
`search collection 'pubchem' for 'Ibuprofen' show (data)` <br> 

- Search for patents which mention a specific smiles molecule: <br> 
`search collection 'patent-uspto' for '\"smiles#ccc(coc(=o)cs)(c(=o)c(=o)cs)c(=o)c(=o)cs\"' show (data)` <br><br>

`display collection matches for '<search_string>' [ save as '<filename.csv>' ]`{: .cmd }
Search all collections for documents that contain a given Deep Search `<search_string>`. This is useful when narrowing down document collection(s) for subsequent search. You can use the `<index_key>` from the returned table in your next search. <br> 

Use the `save as` clause to save the results as a csv file in your current workspace. <br> 

Example: <br> 
`display collection matches for 'Ibuprofen'` <br><br>

<br>

### Collections

`display collections in domains from list <list_of_domains> [ save as '<filename.csv>' ]`{: .cmd }
Display collections that belong to the listed domains. <br> 

Use the `save as` clause to save the results as a csv file in your current workspace. <br> 

Use the command `display all collections` to find available domains. <br> 

Example: <br> 
`display collections in domains from list ['Scientific Literature']` <br><br>

`display all collections [ save as '<filename.csv>' ]`{: .cmd }
Display all available collections in Deep Search. <br> 

Use the `save as` clause to save the results as a csv file in your current workspace. <br><br>

`display collections for domain '<domain_name>'`{: .cmd }
Display the available collections in a given Deep Search domain. <br> 

Use the command `display all collections` to find available domains. <br> 

Example: <br> 
`display collections for domain 'Business Insights'` <br><br>

`display collection details '<collection_name_or_key>'`{: .cmd }
Display the details for a specified collection. You can specify a collection by its name or key. <br> 

Use the command `display all collections` to list available collections. <br> 

Example: <br> 
`display collection details 'Patents from USPTO'` <br><br>

<br>

</details>

## RXN


<details markdown="block">
<summary>See commands</summary>

### General

`interpret recipe '<recipe_paragraph>' | '<txt_filename>'`{: .cmd }
Build a ordered list of actions interpreted from a provided text-based recipe. The recipe can be provided as a string or as a text file from your current workspace. <br> 

Examples: <br> 
- `interpret recipe 'my_recipe.txt'` <br> 
- `interpret recipe 'A solution of ((1S,2S)-1-{[(methoxymethyl-biphenyl-4-yl)-(2-pyridin-2-yl-cyclopropanecarbonyl)-amino]-methyl}-2-methyl-butyl)-carbamic acid tert-butyl ester (25 mg, 0.045 mmol) and dichloromethane (4 mL) was treated with a solution of HCl in dioxane (4 N, 0.5 mL) and the resulting reaction mixture was maintained at room temperature for 12 h. The reaction was then concentrated to dryness to afford (1R,2R)-2-pyridin-2-yl-cyclopropanecarboxylic acid ((2S,3S)-2-amino-3-methylpentyl)-(methoxymethyl-biphenyl-4-yl)-amide (18 mg, 95% yield) as a white solid.'` <br><br>

`list rxn models`{: .cmd }
Lists all RXN AI models currently available. <br><br>

<br>

### Retrosynthesis

`predict retrosynthesis '<smiles>' [ using (option1=<value> option2=<value>) ]`{: .cmd }
Perform a retrosynthesis route prediction on a molecule. <br> 

Optional Parameters that can be specified in the `using` clause: <br> 
- `availability_pricing_threshold=<int>` Maximum price in USD per g/ml of compounds. Default: no threshold. <br> 
- `available_smiles='<smiles>.<smiles>.<smiles>'` List of molecules available as precursors, delimited with a period. <br> 
- `exclude_smiles='<smiles>.<smiles>.<smiles>'` List of molecules to exlude from the set of precursors, delimited with a period. <br> 
- `exclude_substructures='<smiles>.<smiles>.<smiles>'` List of substructures to excludefrom the set of precursors, delimited with a period. <br> 
- `exclude_target_molecule=<boolean>` Excluded target molecule. The default is True <br> 
- `fap=<float>` Every retrosynthetic step is evaluated with the FAP, and is only retained when forward confidence is greater than the FAP value. The default is 0.6. <br> 
- `max_steps=<int>` The maximum number steps in the results. The default is 3. <br> 
- `nbeams=<int>` The maximum number of beams exploring the hypertree. The default is 10. <br> 
- `pruning_steps=<int>` The number of steps to prune a hypertree. The default is 2. <br> 
- `ai_model='<model_name>'` What model to use. Use the command `list rxn models` to list all available models. The default is '2020-07-01'. <br> 

Example: <br> 
`predict retrosynthesis 'BrCCc1cccc2c(Br)c3ccccc3cc12' using (max_steps=3)` <br><br>

<br>

### Prediction

`predict reaction in batch from dataframe <dataframe_name> | file '<filename.csv>' | list ['<smiles>.<smiles>','<smiles>.<smiles>'] [ using (ai_model='<ai_model>') ] [ use_saved ]`{: .cmd }
Run a batch of reaction predictions. The provided list of reactions can be specified as a DataFrame, a CSV file from your current workspace or a list of strings. When proving a DataFrame or CSV file, we will look for the "reactions" column. <br> 

Reactions are defined by combining two SMILES strings delimited by a period. For example: `'BrBr.c1ccc2cc3ccccc3cc2c1'` <br> 

Optional Parameters that can be specified in the `using` clause: <br> 
- `ai_model='<model_name>'` What model to use. Use the command `list rxn models` to list all available models. The default is '2020-07-01'. <br> 

You can reuse previously generated results by appending the optional `use_saved` clause. This will reuse the results of a previously run command with the same parameters, if available. <br> 

Examples: <br> 
- `predict reaction in batch from list ['BrBr.c1ccc2cc3ccccc3cc2c1CCO' , 'BrBr.c1ccc2cc3ccccc3cc2c1']` <br> 
- `predict reaction in batch from list ['BrBr.c1ccc2cc3ccccc3cc2c1CCO' , 'BrBr.c1ccc2cc3ccccc3cc2c1'] use_saved` <br><br>

`predict reaction '<smiles>.<smiles>' [ using (ai_model='<ai_model>') ] [ use_saved ]`{: .cmd }
Predict the reaction between two molecules. <br> 

Reactions are defined by combining two SMILES strings delimited by a period. For example: `'BrBr.c1ccc2cc3ccccc3cc2c1'` <br> 

Optional Parameters that can be specified in the `using` clause: <br> 
- `ai_model='<model_name>'` What model to use. Use the command `list rxn models` to list all available models. The default is '2020-07-01'. <br> 

You can reuse previously generated results by appending the optional `use_saved` clause. This will reuse the results of a previously run command with the same parameters, if available. <br> 

Examples: <br> 
- `predict reaction 'BrBr.c1ccc2cc3ccccc3cc2c1CCO'` <br> 
- `predict reaction 'BrBr.c1ccc2cc3ccccc3cc2c1CCO' use_saved` <br><br>

`predict reaction topn in batch from dataframe <dataframe_name> | file '<filename.csv>' | list ['<smiles>.<smiles>','<smiles>.<smiles>'] [ using (topn=<integer> ai_model='<ai_model>') ] [ use_saved ]`{: .cmd }
Run a batch of reaction predictions for topn. The provided list of reactions can be specified as a DataFrame, a CSV file from your current workspace or a list of strings. When proving a DataFrame or CSV file, we will look for the "reactions" column. <br> 

Reactions are defined by combining two SMILES strings delimited by a period. For example: `'BrBr.c1ccc2cc3ccccc3cc2c1'` <br> 

Optional Parameters that can be specified in the `using` clause: <br> 
- `ai_model='<model_name>'` What model to use. Use the command `list rxn models` to list all available models. The default is '2020-07-01'. <br> 
- `topn=<integer>` Defined the number of results being returned. The default value is 3. <br> 

You can reuse previously generated results by appending the optional `use_saved` clause. This will reuse the results of a previously run command with the same parameters, if available. <br> 

Examples: <br> 
- `predict reaction topn in batch from list ['BrBr.c1ccc2cc3ccccc3cc2c1CCO' , 'BrBr.c1ccc2cc3ccccc3cc2c1']` <br> 
- `predict reaction topn in batch from list ['BrBr.c1ccc2cc3ccccc3cc2c1CCO' , 'BrBr.c1ccc2cc3ccccc3cc2c1'] using (topn=6)` <br> 
- `predict reaction topn in batch from list ['BrBr.c1ccc2cc3ccccc3cc2c1CCO' , 'BrBr.c1ccc2cc3ccccc3cc2c1'] use_saved ` <br><br>

<br>

</details>

## ST4SD


<details markdown="block">
<summary>See commands</summary>

</details>
