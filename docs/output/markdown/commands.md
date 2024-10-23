---
title: Commands
layout: home
nav_order: 4
---

<!--

DO NOT EDIT
-----------
This file auto-generated.
To update it, consult instructions:
https://github.com/acceleratedscience/open-ad-toolkit/tree/main/docs

-->

### Table of Contents
- [OpenAD](#openad)
  - [Small Molecules](#small-molecules)
  - [Molecule Working Set](#molecule-working-set)
  - [Macromolecules](#macromolecules)
  - [General](#general)
  - [Workspaces](#workspaces)
  - [Toolkits](#toolkits)
  - [Runs](#runs)
  - [Utility](#utility)
  - [GUI](#gui)
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

### Small Molecules

`display molecule|mol <name> | <smiles> | <inchi> | <inchikey> | <cid>`{: .cmd }
Display a molecule's details. <br> 

A molecule's details include its identifiers, synonyms, properties and any analysis results it has been enriched with. <br> 

Supported molecule identifiers: <br> 
- `name` / `synonym` <br> 
- `SMILES` <br> 
- `InChI` <br> 
- `InChIKey` <br> 
- `PubChem CID` <br> 

Notes: <br> 
- You can use the 'mol' shorthand instead of 'molecule'. <br> 
- If the requested molecule exists in your current working set or in memory, that version will be prioritized. <br> 

Examples: <br> 
- Display a molecule by name: <br> 
  `display molecule Aspirin` <br> 

- Display a molecule by SMILES string: <br> 
  `display molecule CC(=O)OC1=CC=CC=C1C(=O)O` <br> 

- Display a molecule by InChI string: <br> 
  `display mol InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)` <br> 

- Display a molecule by InChIKey: <br> 
  `display mol BSYNRYMUTXBXSQ-UHFFFAOYSA-N` <br> 

- Display a molecule by its PubChem CID: <br> 
  `display mol 2244` <br><br>

`show molecule|mol <name> | <smiles> | <inchi> | <inchikey> | <cid>`{: .cmd }
Launch the molecule viewer to visualize your molecule and inspect its properties. <br> 

Supported molecule identifiers: <br> 
- `name` / `synonym` <br> 
- `SMILES` <br> 
- `InChI` <br> 
- `InChIKey` <br> 
- `PubChem CID` <br> 

Notes: <br> 
- You can use the 'mol' shorthand instead of 'molecule'. <br> 
- If the requested molecule exists in your current working set or in memory, that version will be prioritized. <br> 

Examples: <br> 
- Inspect a molecule by name: <br> 
  `show mol aspirin` <br> 

- Inspect a molecule by SMILES string: <br> 
  `show mol CC(=O)OC1=CC=CC=C1C(=O)O` <br> 

- Inspect a molecule by InChI string: <br> 
  `show mol InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)` <br> 

- Inspect a molecule by PubChem CID: <br> 
  `show mol 2244` <br><br>

`show molset|molecule-set '<molset_or_sdf_or_smi_path>' | using dataframe <dataframe>`{: .cmd }
Launch the molset viewer to visualize the contents of a molecule set. <br> 

Examples: <br> 
- `show molset 'neurotransmitters.mol.json'` <br> 
- `show molset 'neurotransmitters.sdf'` <br> 
- `show molset 'neurotransmitters.smi'` <br> 
- `show molset my_dataframe` <br><br>

`@(<name> | <smiles> | <inchi> | <inchikey> | <cid>)>><molecule_property_name>`{: .cmd }
Request a molecule's property. <br> 

Supported molecule identifiers: <br> 
- `name` / `synonym` <br> 
- `SMILES` <br> 
- `InChI` <br> 
- `InChIKey` <br> 
- `PubChem CID` <br> 

Notes: <br> 
- In addition to a molecule's properties or identifiers, you can also request its synonyms. <br> 
- If the requested molecule exists in your current working set or in memory, that version will be prioritized. <br> 
- In addition to the properties listed below, you can also request any additional properties that are available in your molecule working set. <br> 

Examples: <br> 
- Obtain the molecular weight of a molecule known as Aspirin. <br> 
  `@aspirin>>molecular_weight` <br> 

- Obtain the canonical smiles string for the dopamine neurotransmitter. <br> 
  `@dopamine>>canonical_smiles` <br> 

- Obtain a molecules xlogp value using a SMILES string. <br> 
  `@CC(=O)OC1=CC=CC=C1C(=O)O>>xlogp` <br> 

- Obtain all the synonyms for ibuprofen. <br> 
 `@ibuprofen>>synonyms` <br> 

Available properties that can be queried: <br> 
atom_stereo_count                       bond_stereo_count                       canonical_smiles                         <br> 
charge                                  cid                                     complexity                               <br> 
conformer_count_3d                      conformer_id_3d                         conformer_model_rmsd_3d                  <br> 
conformer_rmsd_3d                       coordinate_type                         covalent_unit_count                      <br> 
defined_atom_stereo_count               defined_bond_stereo_count               effective_rotor_count_3d                 <br> 
exact_mass                              feature_acceptor_count_3d               feature_anion_count_3d                   <br> 
feature_cation_count_3d                 feature_count_3d                        feature_donor_count_3d                   <br> 
feature_hydrophobe_count_3d             feature_ring_count_3d                   h_bond_acceptor_count                    <br> 
h_bond_donor_count                      heavy_atom_count                        inchi                                    <br> 
inchikey                                isomeric_smiles                         isotope_atom_count                       <br> 
iupac_name                              mmff94_energy_3d                        mmff94_partial_charges_3d                <br> 
molecular_formula                       molecular_weight                        molecular_weight_exact                   <br> 
monoisotopic_mass                       multipoles_3d                           multipoles_3d                            <br> 
pharmacophore_features_3d               pharmacophore_features_3d               rotatable_bond_count                     <br> 
sol                                     sol_classification                      tpsa                                     <br> 
undefined_atom_stereo_count             undefined_bond_stereo_count             volume_3d                                <br> 
x_steric_quadrupole_3d                  xlogp                                   y_steric_quadrupole_3d                   <br><br>

`export molecule|mol <name> | <smiles> | <inchi> | <inchikey> | <cid> [ as file ]`{: .cmd }
Either save a molecule to your workspace as a JSON file (CLI) or return a dictionary of the molecule's properties (Jupyter Notebook). <br> 

Supported molecule identifiers: <br> 
- `name` / `synonym` <br> 
- `SMILES` <br> 
- `InChI` <br> 
- `InChIKey` <br> 
- `PubChem CID` <br> 

Notes: <br> 
- The requested molecule does not need to be in your current working set. <br> 
- The `as file` clause is only needed when you wish to save a file to your workspace from within a Jupytyer Notebook. <br> 
- You can use the 'mol' shorthand instead of 'molecule'. <br> 
- If the requested molecule exists in your current working set or in memory, that version will be prioritized. <br> 

Examples: <br> 
- Export a molecule by name: <br> 
  `export molecule aspirin` <br> 

- Export a molecule by SMILES string, and force it to save as a file (only relevant in a Jupyter Notebook): <br> 
  `export mol aspirin as file` <br><br>

<br>

### Molecule Working Set

`add molecule|mol <name> | <smiles> | <inchi> | <inchikey> | <cid> [ as '<name>' ] [ basic ] [ force ]`{: .cmd }
Add a molecule to your current molecule working set. <br> 

Supported molecule identifiers: <br> 
- `name` / `synonym` <br> 
- `SMILES` <span style="color: #ccc">- supports [ basic ]</span> <br> 
- `InChI` <span style="color: #ccc">- supports [ basic ]</span> <br> 
- `InChIKey` <br> 
- `PubChem CID` <br> 

Options: <br> 
- `as <name>`: Provide a custom name for the molecule, which will be used by the software whenever refering to it going forward. <br> 
  Note: you can always update a molecule's name later by running `rename molecule <name>`. <br> 
- `basic`: Create a minimal molecule without enriching it with PubChem data. This is only relevant when using a SMILES or InChI string as identifier. Because no API calls are made, this is much faster than the default behavior. <br> 
- `force`: Suppress the confirmation step before adding a molecule, which may be desired in batch operations. <br> 

Notes: <br> 
- You can use the 'mol' shorthand instead of 'molecule'. <br> 
- If the requested molecule exists in your current working set or in memory, that version will be prioritized. <br> 

Examples: <br> 
- Add a molecule by SMILES string: <br> 
  `add molecule CC(=O)OC1=CC=CC=C1C(=O)O` <br> 

- Add a molecule by SMILES string, without enriching it with PubChem data: <br> 
  `add molecule CC(=O)OC1=CC=CC=C1C(=O)O basic` <br> 

- Add a molecule by SMILES string, giving it a custom name: <br> 
  `add molecule CC(=O)OC1=CC=CC=C1C(=O)O as 'mymol' basic` <br> 

- Add a molecule by unquoted InChI string: <br> 
  `add mol InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)` <br> 

- Add a molecule by quoted InChI string: <br> 
  `add mol 'InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)'` <br> 

- Add a molecule by InChIKey: <br> 
  `add mol BSYNRYMUTXBXSQ-UHFFFAOYSA-N` <br> 

- Add a molecule by InChIKey, giving it a custom name: <br> 
  `add mol BSYNRYMUTXBXSQ-UHFFFAOYSA-N as 'mymol'` <br> 

- Add a molecule by name: <br> 
  `add molecule aspirin` <br> 

- Add a molecule by name, supressing the confirmation step: <br> 
  `add molecule aspirin force` <br> 

- Add a molecule by quoted name: <br> 
  `add mol 'Aspirin 325 mg'` <br> 

- Add a molecule by its PubChem CID: <br> 
  `add mol 2244` <br><br>

`remove molecule|mol <name> | <smiles> | <inchi> | <inchikey> | <cid> [ force ]`{: .cmd }
Remove a molecule from the current working set based on a given molecule identifier. <br> 

Supported molecule identifiers: <br> 
- `name` / `synonym` <br> 
- `SMILES` <br> 
- `InChI` <br> 
- `InChIKey` <br> 
- `PubChem CID` <br> 

Options: <br> 
- `force`: Suppress the confirmation step before removing a molecule, which may be desired in batch operations. <br> 

Notes: <br> 
- You can use the 'mol' shorthand instead of 'molecule'. <br> 

Examples: <br> 
- Remove a molecule by name: <br> 
  `remove molecule Aspirin` <br> 

- Remove a molecule by SMILES: <br> 
  `remove molecule CC(=O)OC1=CC=CC=C1C(=O)O` <br> 

- Remove a molecule by InChIKey: <br> 
  `remove mol BSYNRYMUTXBXSQ-UHFFFAOYSA-N` <br> 

- Remove a molecule by InChI <br> 
  `remove mol  InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)` <br> 

- Remove a molecule by CID: <br> 
  `remove mol 2244` <br><br>

`list molecules|mols`{: .cmd }
List all molecules in the current working set. <br> 

Notes: <br> 
- You can use the 'mols' shorthand instead of 'molecules'. <br><br>

`show molecules|mols`{: .cmd }
Launch the molset viewer to visualize your molecule working set. <br> 

Notes: <br> 
- You can use the 'mols' shorthand instead of 'molecules'. <br><br>

`enrich molecules|mols with analysis`{: .cmd }
Enrich the molecules in your current working set with the results of the last performed analysis. <br> 

This assumes that your molecule working set contains either the input molecule or any of the result molecules from the analysis. <br> 

Notes: <br> 
- You can use the 'mols' shorthand instead of 'molecules'. <br> 

Currently supported analysis commands: <br> 

RXN: <br> 
- `predict reaction` <br> 
- `predict retrosynthesis` <br> 

DS4SD: <br> 
- `search for patents containing molecule` <br> 
- `search for similiar molecules` <br> 

Please refer to the DS4SD and RXN toolkits for further assistance on these commands. For example: <br> 
`set context rxn` <br> 
`predict reaction ?` <br><br>

`clear analysis cache`{: .cmd }
Clear the analysis results cache for your current workspace. <br> 

Please refer to the `enrich molecules|mols with analysis` command for more information about analysis results. <br><br>

`display sources <name> | <smiles> | <inchi> | <inchikey> | <cid>`{: .cmd }
Display the sources of a molecule's properties, attributing how they were calculated or where they were sourced. <br> 

Supported molecule identifiers: <br> 
- `name` / `synonym` <br> 
- `SMILES` <br> 
- `InChI` <br> 
- `InChIKey` <br> 
- `PubChem CID` <br> 

Notes: <br> 
- If the requested molecule exists in your current working set or in memory, that version will be prioritized. <br><br>

`rename molecule|mol <molecule_identifer_string> as <molecule_name>`{: .cmd }
Rename a molecule in the current working set. <br> 

Notes: <br> 
- You can use the 'mol' shorthand instead of 'molecule'. <br> 

Example: <br> 
- Assuming you've added the molecule <span style="color: #dc0">CC(=O)OC1=CC=CC=C1C(=O)O</span> to your molecule working set, you can then rename it as such: <br> 
  `rename molecule CC(=O)OC1=CC=CC=C1C(=O)O as Aspirin` <br><br>

`load molecules|mols from file '<filename.molset.json|sdf|csv|smi>' [ enrich ] [ append ]`{: .cmd }
Load molecules from a file into your molecule working set. <br> 

Supported file formats: <br> 
- molset (.molset.json) <br> 
- SDF (.sdf) <br> 
- CSV (.csv) <br> 
- SMILES (.smi) <br> 

Options: <br> 
- Append `enrich` to enrich the molecule with data from pubchem. <br> 
- Append `append` to append the molecules to the existing working set instead of overwriting it. <br> 

Notes: <br> 
- To see some example input files, you can export the molecules from your working set using the `export molecules|mols as ...` command. <br> 
- You can use the 'mols' shorthand instead of 'molecules'. <br> 

Examples: <br> 
- Load molecules from a molset JSON file: <br> 
  `load molecules from file 'my_molecules.molset.json'` <br> 

- Load molecules from an SDF file, appending them to the existing working set: <br> 
  `load mols from file 'my_molecules.sdf'` append` <br> 

- Load molecules from an CSV file: <br> 
  `load molecules from file 'my_molecules.csv'` <br> 

- Load molecules from an SMILES file, enriching them with PubChem data: <br> 
  `load mols from file 'my_molecules.smi'` enrich` <br><br>

`load molecules|mols from dataframe <dataframe> [ enrich ] [ append ]`{: .cmd }
Load molecules from a dataframe into your molecule working set. <br> 

Options: <br> 
- Append `enrich` to enrich the molecule with data from pubchem. <br> 
- Append `append` to append the molecules to the existing working set instead of overwriting it. <br> 

Notes: <br> 
- This command only works when called from a Jupyter Notebook or the API. <br> 
- To see some example input files, you can export the molecules from your working set using the `export molecules|mols as ...` command. <br> 
- You can use the 'mols' shorthand instead of 'molecules'. <br> 

Examples: <br> 
- Load molecules from a dataframe: <br> 
  `load molecules from dataframe my_dataframe` <br> 

- Load molecules from a dataframe, appending them to the existing working set: <br> 
  `load mols from dataframe my_dataframe append` <br> 

- Load molecules from a dataframe, enriching them with PubChem data: <br> 
  `load mols from dataframe my_dataframe enrich` <br><br>

`merge molecules|mols data from dataframe <dataframe> [ enrich ]`{: .cmd }
Merge molecule data from a dataframe into the molecules in your working set. <br> 

Options: <br> 
- Append `enrich` to enrich the molecule with data from pubchem. <br> 

The dataframe columns should be named as follows: <br> 
- `subject` or `smiles`: molecules similes string <br> 
- `property`: the name of the property to be merged <br> 
- `result`: the value of the property to be nmerged <br> 

Sample input file: <br> 

 subject                                                               property                        result  <br> 
 --------------------------------------------------------------------  -------------------------  -----------  <br> 
 O=C(N)C(F)(OC(F)(F)C(F)(F)C(F)(F)F)C(F)(F)F                           molecular_weight               329.065  <br> 
 O=C(O)C(F)(NC(F)(F)C(F)(F)C(F)(F)F)C(F)(F)F                           molecular_weight               329.065  <br> 
 O=C(O)C(F)(OC(F)(F)C(F)(F)C)CF                                        molecular_weight               240.099  <br> 
 O=C(O)C(F)(OC(O)(F)C(F)(F)C(F)(F)F)C(F)(F)F                           molecular_weight               328.058  <br> 
 O=C(O)C(F)(OC(Cl)(F)C(F)(F)C(F)(F)F)C(F)(F)F                          molecular_weight               346.504  <br> 
 O=C(O)C(F)(OC(F)(F)C(F)(O)C(F)(F)F)C(F)(F)F                           molecular_weight               328.058  <br> 
 O=C(O)C(F)(OC(F)(O)C(F)(F)C(F)(F)F)C(F)(F)F                           molecular_weight               328.058  <br> 
 O=C(O)C(F)(OC(F)(F)C(F)(Br)C(F)(F)F)C(F)(F)F                          molecular_weight               390.955  <br> 
 O=C(O)C(F)OC(O)(F)C(F)(F)C(F)(F)F                                     molecular_weight               260.061  <br> 

Notes: <br> 
- You can use the 'mol' shorthand instead of 'molecule'. <br> 

Examples: <br> 
- Merge molecule data from a dataframe called `new_props`: <br> 
  `merge molecules data from dataframe new_props` <br> 

- Merge molecule data from a dataframe called `new_props`, while enriching the molecules with PubChem data: <br> 
  `merge molecules data from dataframe new_props enrich` <br><br>

`export molecules|mols [ as '<filename.molset.json|sdf|csv|smi>' ]`{: .cmd }
Export your molecule working set as a file (CLI) or return it as a dataframe (Jupyter/API). <br> 

Supported file formats: <br> 
- molset (.molset.json) <br> 
- SDF (.sdf) <br> 
- CSV (.csv) <br> 
- SMILES (.smi) <br> 

Notes: <br> 
- When exporting as a file, the filename's extension will define what format the molecule are exported as. <br> 
- If no filename or extension is provided, the molecules will be saved as CSV file. <br> 
- When run inside a Jupyter Notebook or from the API, the `as <filename>` clause will be ignored and a dataframe will be returned. <br> 
- You can use the 'mols' shorthand instead of 'molecules'. <br><br>

`clear molecules|mols [ force ]`{: .cmd }
Clear the molecule working set. <br> 

Options: <br> 
- `force`: Suppress the confirmation step before clearing the working set, which may be desired in batch operations. <br> 

Notes: <br> 
- You can use the 'mols' shorthand instead of 'molecules'. <br><br>

<br>

### Macromolecules

`show mmol|protein <fasta> | '<pdb_id>'`{: .cmd }
Launch the molecule viewer to visualize your macromolecule and inspect its properties. <br> 

Examples: <br> 
- Show a protein by its PDBe ID: <br> 
  `show mmol '2g64'` <br> 

- Show a protein by its FASTA string: <br> 
  `show protein MAKWVCKICGYIYDEDAGDPDNGISPGTKFEELPDDWVCPICGAPKSEFEKLED` <br><br>

<br>

### General

`openad`{: .cmd }
Display the openad splash screen. <br><br>

`get status`{: .cmd }
Display the currently selected workspace and toolkit. <br><br>

`display history`{: .cmd }
Display the last 30 commands run in your current workspace. <br><br>

`clear sessions [ force ]`{: .cmd }
Clear any other sessions that may be running. <br><br>

<br>

### Workspaces

`list workspaces`{: .cmd }
Lists all your workspaces. <br><br>

`set workspace <workspace_name>`{: .cmd }
Change the current workspace. <br><br>

`get workspace [ <workspace_name> ]`{: .cmd }
Display details a workspace. When no workspace name is defined, details of your current workspace are displayed. <br><br>

`get workspace path`{: .cmd }
Returns your current workspace path. This is useful when you wish to access your workspace files from a Jupyter Notebook. <br><br>

`show workspace`{: .cmd }
Open the file browsing interface displaying your workspace files and directories. <br><br>

`open workspace [ <workspace_name> ]`{: .cmd }
Open your workspace directory in the file explorer of your operating system. This is useful when you want to move or rename files. When no workspace name is defined, your current workspace will open. <br><br>

`create workspace <workspace_name> [ description('<description>') on path '<path>' ]`{: .cmd }
Create a new workspace with an optional description and path. <br><br>

`remove workspace <workspace_name>`{: .cmd }
Remove a workspace from your registry. Note that this doesn't remove the workspace's directory. <br><br>

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

`remove run <run_name>`{: .cmd }
remove a run. <br><br>

`save run as <run_name>`{: .cmd }
Stop recording a run and save it. <br><br>

`run <run_name>`{: .cmd }
Execute a previously recorded run. This will execute every command and continue regardless of any failures. <br><br>

`list runs`{: .cmd }
List all runs saved in the current workspace. <br><br>

`display run <run_name>`{: .cmd }
Display the commands stored in a certain run. <br><br>

<br>

### Utility

`display data '<filename.csv>'`{: .cmd }
Display data from a csv file. <br><br>

`-> result save [as '<filename.csv>']`{: .cmd }
Save table data to csv file. <br><br>

`-> result open`{: .cmd }
Explore table data in the browser. <br> 
        if you append `-d` to the end of the command `result open -d` display will result to data viewer. <br><br>

`-> result edit`{: .cmd }
Edit table data in the browser. <br> 
        if you append `-d` to the end of the command `result open -d` display will result to data viewer. <br><br>

`-> result copy`{: .cmd }
Copy table data to clipboard, formatted for spreadheet. <br><br>

`-> result display`{: .cmd }
Display the result in the CLI. <br> 

        if you append `-d` to the end of the command `result open -d` display will result to data viewer. <br><br>

`-> result as dataframe`{: .cmd }
Return the result as dataframe (only for Jupyter Notebook) <br><br>

`edit config '<json_config_file>' [ schema '<schema_file>']`{: .cmd }
Edit any JSON file in your workspace directly from the CLI. If a schema is specified, it will be used for validation and documentation. <br><br>

<br>

### GUI

`install gui`{: .cmd }
Install the OpenAD GUI (graphical user interface). <br> 

The graphical user interface allows you to browse your workspace and visualize your datasets and molecules. <br><br>

`launch gui`{: .cmd }
Launch the OpenAD GUI (graphical user interface). <br><br>

`restart gui`{: .cmd }
Terminate and then restart the GUI server. <br><br>

`quit gui`{: .cmd }
Terminate the GUI server. <br><br>

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

`list files [ path ]`{: .cmd }
List al directories and files in your current workspace. <br><br>

`import from '<external_source_file>' to '<workspace_file>'`{: .cmd }
Import a file from outside OpenAD into your current workspace. <br><br>

`export from '<workspace_file>' to '<external_file>'`{: .cmd }
Export a file from your current workspace to anywhere on your hard drive. <br><br>

`copy file '<workspace_file>' to '<other_workspace_name>'`{: .cmd }
Export a file from your current workspace to another workspace. <br><br>

`remove '<filename>'`{: .cmd }
Remove a file from your current workspace. <br><br>

`open '<filename>'`{: .cmd }
Open a file or dataframe in an iframe  <br> 

Examples: <br> 
- `open 'base_molecules.sdf'` <br> 
- `open my_dataframe` <br><br>

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

`model auth add group '<auth_group>'|<auth_group> with '<api_key>'`{: .cmd }
add an authentication group for model services to use <br><br>

`model auth remove group '<auth_group>' | <auth_group>`{: .cmd }
remove an authentication group <br><br>

`model auth add service '<service_name>'|,service_name> to group '<auth_group>'|<auth_group>`{: .cmd }
Attach an authentication group to a model service <br><br>

`model auth remove service '<service_name>'|<service_name>`{: .cmd }
Detatch an authentication group from a model service <br><br>

`model service status`{: .cmd }
Get the status of currently cataloged services <br><br>

`model service describe '<service_name>'|<service_name>`{: .cmd }
get the configuration of a service <br><br>

`model catalog list`{: .cmd }
get the list of currently cataloged services <br><br>

`uncatalog model service '<service_name>'|<service_name>`{: .cmd }
uncatalog a model service  <br> 

 Example:  <br> 
`uncatalog model service 'gen'` <br><br>

`catalog model service from (remote) '<path> or <github> or <service_url>' as  '<service_name>'|<service_name>   USING (<parameter>=<value> <parameter>=<value>)`{: .cmd }
catalog a model service from a path or github or remotely from an existing OpenAD service. <br> 
(USING) optional headers parameters for communication with service backend. <br> 
If you are cataloging a service using a model defined in a directory, provide the absolute ` <path> ` of that directory in quotes. <br> 

The following options require the `remote` option be declared. <br> 

If you are cataloging a service using a model defined in github repository, provide the absolute ` <github> ` of that github directory quotes. <br> 

If you are cataloging a remote service on a ip address and port provide the remote services ipaddress and port in quoted string e.g. `'0.0.0.0:8080'` <br> 

`service_name`: this is the name of the service as you will define it for your usage. e.g `prop` short for properties.  <br> 

USING Parameters: <br> 

If using a hosted service the following parameters must be supplied: <br> 
-`Inference-Service`: this is the name of the inference service that is hosted, it is a required parameter if cataloging a remote service. <br> 
An authorization parameter is always required if cataloging a hosted service, either Auhtorisation group (`auth_group`) or Authorisation bearer_token/api_key (`Authorization`): <br> 
-`auth_group`: this is the name of an authorization group which contains the api_key linked to the service access. This can only be used if `Authorization` is not also defined. <br> 
OR <br> 
-`Authorization`: this parameter is designed to be used when a `auth_group` is not defined. <br> 

Example: <br> 

Skypilot Deployment <br> 
-`catalog model service from 'git@github.com:acceleratedscience/generation_inference_service.git' as 'gen'` <br> 

Service using a authentication group  <br> 
-`catalog model service from remote '<service_url>' as  molf  USING (Inference-Service=molformer  )` <br> 
` model auth add service 'molf' to group 'default'` <br> 

Single Authorisation Service <br> 
-`openad catalog model service from remote '<service_URL>' as 'gen' USING (Inference-Service=generation Authorization='<api_key>')` <br> 

Catalog a remote service shared with you: <br> 
-`catalog model service from remote 'http://54.235.3.243:30001' as gen` <br><br>

`model service up '<service_name>'|<service_name> [no_gpu]}`{: .cmd }
launches a cataloged model service when it was cataloged as a self managed service from a directory or github repository. <br> 
If you do not want to launch a service with GPU you should specify `no_gpu` at the end of the command. <br> 
Examples: <br> 

-`model service up gen` <br> 

-`model service up 'gen'` <br> 

-`model service up gen no_gpu` <br><br>

`model service local up '<service_name>'|<service_name>`{: .cmd }
Launches a model service locally. <br> 

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

RXN was trained on more than 3 million chemical reactions, derived from publicly available patents. Since then, the Molecular Transformer has outperformed all data-driven models, achieving more than 90% accuracy on forward chemical reaction predictions (reactants + reagents to products) <br> 

  > **Note:** The `using` clause requires all enclosed parameters to be defined in the same order as listed below. <br> 

Optional Parameters that can be specified in the `using` clause: <br> 
- `availability_pricing_threshold=<int>` Maximum price in USD per g/ml of compounds. Default: no threshold. <br> 
- `available_smiles='<smiles>.<smiles>.<smiles>'` List of molecules available as precursors, delimited with a period. <br> 
- `exclude_smiles='<smiles>.<smiles>.<smiles>'` List of molecules to exclude from the set of precursors, delimited with a period. <br> 
- `exclude_substructures='<smiles>.<smiles>.<smiles>'` List of substructures to exclude from the set of precursors, delimited with a period. <br> 
- `exclude_target_molecule=<boolean>` Excluded target molecule. The default is True <br> 
- `fap=<float>` Every retrosynthetic step is evaluated with the FAP, and is only retained when forward confidence is greater than the FAP value. The default is 0.6. <br> 
- `max_steps=<int>` The maximum number steps in the results. The default is 3. <br> 
- `nbeams=<int>` The maximum number of beams exploring the hypertree. The default is 10. <br> 
- `pruning_steps=<int>` The number of steps to prune a hypertree. The default is 2. <br> 
- `ai_model='<model_name>'` What model to use. Use the command `list rxn models` to list all available models. The default is '2020-07-01'. <br> 

There are different models available for use with this command including: '12class-tokens-2021-05-14', '2019-09-12', '2020-04-24', '2020-07-01', '2020-07-31', 'aizynth-2020-06-11', 'disconnection-aware-2022-06-24', 'enzymatic-2021-04-16', 'enzymatic-2022-05-31', 'sulfonium-2020-10-27' <br> 

Examples: <br> 
`predict retrosynthesis 'BrCCc1cccc2c(Br)c3ccccc3cc12' using (max_steps=3)` <br> 

`predict retrosynthesis  'BrCCc1cccc2c(Br)c3ccccc3cc12' using (max_steps=6 ai_model='12class-tokens-2021-05-14' ) ` <br><br>

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

