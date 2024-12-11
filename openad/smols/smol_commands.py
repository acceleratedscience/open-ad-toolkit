"""Functions that are called for molecule commands"""

import os
import json
import glob
import copy
import pickle
import urllib.parse
from copy import deepcopy
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

# Helpers
from openad.helpers.output import output_text, output_table, output_warning, output_error, output_success
from openad.helpers.output_msgs import msg
from openad.helpers.pretty_data import list_columns, key_val_columns, key_val_full
from openad.app.global_var_lib import GLOBAL_SETTINGS


# Molecule functions
from openad.smols.smol_transformers import molset2dataframe
from openad.smols.smol_functions import (
    SMOL_PROPERTIES,
    PCY_IDFR,
    find_smol,
    get_smol_from_mws,
    get_smol_from_pubchem,
    get_human_properties,
    get_smol_from_list,
    mws_add,
    mws_remove,
    merge_molecule_REPLACE,
    save_molset_as_json,
    save_molset_as_sdf,
    save_molset_as_csv,
    save_molset_as_smiles,
    clear_mws,
)


def display_molecule(cmd_pointer, inp):
    """
    Display a molecule's identifiers, synonyms and properties.
    """

    molecule_identifier = inp.as_dict()["molecule_identifier"]

    smol = find_smol(cmd_pointer, molecule_identifier, show_spinner=True)
    if smol is None:
        output_error("Error: Not a valid Molecule", return_val=False)
        return None

    if smol is not None:
        output = format_identifiers(smol) + "\n\n" + format_synonyms(smol) + "\n\n" + format_properties(smol)
        if "analysis" in smol and smol["analysis"] != []:
            output = output + "\n\n" + format_analysis(smol)

    if GLOBAL_SETTINGS["display"] == "notebook":
        import py3Dmol
        from IPython.display import Markdown, display, HTML

        astyle = "stick"
        view_mol = Chem.MolFromSmiles(smol["identifiers"]["canonical_smiles"])  # pylint: disable=no-member
        view_mol = Chem.AddHs(view_mol)  # pylint: disable=no-member
        AllChem.EmbedMolecule(view_mol)  # pylint: disable=no-member
        AllChem.MMFFOptimizeMolecule(view_mol, maxIters=200)  # pylint: disable=no-member
        mblock = Chem.MolToMolBlock(view_mol)  # pylint: disable=no-member

        view = py3Dmol.view(width=700, height=500)
        view.addModel(mblock, "mol")
        view.setStyle({astyle: {"model": -1}})
        view.zoomTo()
        view.animate({"loop": "forward"})
        view.show()
        output = output.replace("<soft>", "<span style=color:#cccccc;white-space=pre>")
        output = output.replace("</soft>", "</span>")
        output = output.replace("<cyan>", "<span style=color:#00AAAA;white-space=pre>")
        output = output.replace("</cyan>", "</span>")
        output = output.replace("\n", "<br>")
        display(HTML("<pre>" + output + "</pre>"))

        # display(output_text(print_string, return_val=True))
    else:
        output_text(output, edge=True, pad=1)

    return True


def format_identifiers(smol):
    """
    Format the identifiers for display.
    """

    output = f"<h1>{(smol['identifiers'].get('name') or 'Unknown molecule' ).capitalize()}</h1>\n"
    identifiers = smol.get("identifiers", {})
    output = output + key_val_full(identifiers)
    return output


def format_synonyms(smol, truncate=30):
    """
    Format synonyms for display.
    """

    synonyms = smol.get("synonyms", [])
    synonym_count = len(synonyms)

    # Title
    output = "<h1>Synonyms:</h1>\n"

    # No synonyms found
    if synonym_count == 0:
        return output + "\n<soft>No synonyms found.</soft>"

    # Truncate list of synonyms
    if truncate is not None:
        synonyms = synonyms[:truncate]
        is_truncated = synonym_count > len(synonyms)
    else:
        is_truncated = False

    # Truncation note
    if is_truncated:
        output = (
            output
            + f"<soft>This list is truncated. To view all synonyms, run <cmd>@{smol['identifiers']['name']}>>synonyms</cmd></soft>\n\n"
        )

    # Values
    output = output + list_columns(synonyms, is_truncated=is_truncated)

    return output


def format_properties(smol):
    """
    Format properties for display.
    """

    properties = get_human_properties(smol)

    # Title
    output = "<h1>Properties:</h1>\n"

    # No properties found
    if len(properties) == 0:
        return output + "\n<soft>No properties found.</soft>"

    output = output + key_val_columns(properties, ignore_keys=["DS_URL"])
    return output


def format_analysis(smol):
    """
    Format analysis for display.
    """

    if "analysis" not in smol or smol["analysis"] == []:
        return ""

    output = "<h1>Analysis:</h1>"

    for item in smol["analysis"]:
        results = copy.deepcopy(item["results"])
        del item["results"]

        # Compile analysis context header
        header_output = []
        for key, val in item.items():
            val = "<soft>-</soft>" if not val and val != 0 else val
            header_output.append(f"<cyan>{key.capitalize()}:</cyan> {val}")
        header_output = " <soft>/</soft> ".join(header_output)
        output = output + "\n" + header_output

        # TEMPORARY FIX
        # RXN's 'predict retrosynthesis' results are formatted
        # as a numbered object instead of a list. This causes
        # key_val_full and thus display_molecule to fail.
        # Fixing this can happen when RXN is moved to a plugin.
        # Until then, this fix is a workaround.
        # So results formatted as:
        #   { 0: {}, 1: {} }
        # are turned into:
        #   [ {}, {} ]
        if isinstance(results, dict) and int(list(results.keys())[0]) == 0 and int(list(results.keys())[1]) == 1:
            results = [result for result in results.values()]

        # Compile results
        results_output = ""
        for i, result in enumerate(results):
            results_output = results_output + f"\n\n<soft>Result #{i}</soft>" + key_val_full(result, indent=2)

        # Compile final results output
        output = output + results_output
    return output


def display_property_sources(cmd_pointer, inp):
    """displays a molecule properties sources"""

    molecule_identifier = inp.as_dict()["molecule_identifier"]

    mol = get_smol_from_mws(cmd_pointer, molecule_identifier)

    if mol is not None:
        print_string = format_sources(mol)

    else:
        mol = get_smol_from_pubchem(molecule_identifier)
        if mol is not None:
            cmd_pointer.last_external_molecule = mol.copy()
            print_string = format_sources(mol)

        else:
            output_error("Molecule not available on pubchem", return_val=False)
            return None
    if GLOBAL_SETTINGS["display"] == "notebook":
        from IPython.display import Markdown, display, HTML

        print_string = print_string.replace("<success>", "<text style=color:green;white-space=pre>")
        print_string = print_string.replace("</success>", "</text>")
        print_string = print_string.replace("<yellow>", "<b>")
        print_string = print_string.replace("</yellow>", "</b>")
        print_string = print_string.replace("\n", "<br>")
        display(HTML("<pre>" + print_string + "</pre>"))
    else:
        output_text(print_string, edge=True)

    return True


def format_sources(smol):
    """
    Format the sources for display.
    """

    output = f"\n<yellow>Name:</yellow> {smol['identifiers']['name']}\n"

    for mol_property, source in smol["property_sources"].items():
        if mol_property not in smol["properties"] or smol["properties"][mol_property] is None:
            continue
        output = output + "\n\n<yellow>Property:</yellow> {} \n".format(mol_property)
        output = output + key_val_columns(source)

    return output


def export_molecule(cmd_pointer, inp):
    """
    Export a molecule as a smol.json file.
    """

    molecule_identifier = inp.as_dict()["molecule_identifier"]
    smol = find_smol(cmd_pointer, molecule_identifier)
    if smol:
        cmd_pointer.last_external_molecule = smol

    if "as_file" in inp.as_dict() or GLOBAL_SETTINGS["display"] in ["terminal"]:
        json_file = open(
            cmd_pointer.workspace_path(cmd_pointer.settings["workspace"].upper())
            + "/"
            + smol["identifiers"]["name"]
            + ".smol.json",
            "w",
            encoding="utf-8",
        )
        json.dump(smol, json_file)
        output_success(
            "Molecule saved to your workspace as <yellow>" + smol["identifiers"]["name"] + ".smol.json</yellow>.",
            return_val=False,
        )
    elif GLOBAL_SETTINGS["display"] in ["api", "notebook"]:
        return deepcopy(smol)
    return True


def add_molecule(cmd_pointer, inp):
    """
    Adds a molecule to your molecule working set.
    """

    identifier = inp.as_dict()["molecule_identifier"]
    basic = "basic" in inp.as_dict()
    name = inp.as_dict()["name"] if "name" in inp.as_dict() else identifier
    force = "force" in inp.as_dict()

    # Create molecule dict.
    smol = find_smol(cmd_pointer, identifier, name, basic)

    # Add it to the working set.
    if smol:
        mws_add(cmd_pointer, smol, force=force)


def remove_molecule(cmd_pointer, inp):
    """
    Removes a molecule from your molecule working set.
    """

    molecule_identifier = inp.as_dict()["molecule_identifier"]
    force = "force" in inp.as_dict()

    # Remove molecule from working set.
    mol = get_smol_from_mws(cmd_pointer, molecule_identifier)
    if not mol:
        return output_error(f"Molecule <yellow>{molecule_identifier}</yellow> was not found in your working set.")
    mws_remove(cmd_pointer, mol, force=force)


def clear_molecules(cmd_pointer, inp):
    """
    Clear the molecule working set.
    """

    force = "force" in inp.as_dict()
    clear_mws(cmd_pointer, force)


def list_molecules(cmd_pointer, inp):
    """lists all molecules in the working set"""
    display_list = pd.DataFrame()

    if len(cmd_pointer.molecule_list) > 0:
        for smol in cmd_pointer.molecule_list:
            identifiers = smol.get("identifiers", {})
            display_list = pd.concat([display_list, pd.DataFrame([identifiers])])
        return display_list
        # if GLOBAL_SETTINGS["display"] == "notebook":
        #    return output_table(display_list, is_data=True)
        # else:
        #    return output_table(display_list)

    else:
        return output_warning("Your molecule working set is empty")


def show_molecules(cmd_pointer, inp):
    """Display your working set in the browser or iframe"""
    from openad.gui.gui_launcher import gui_init

    gui_init(cmd_pointer, "my-mols")


def rename_mol_in_list(cmd_pointer, inp):
    """
    Rename a molecule in your molecule working set.
    """

    identifier = inp.as_dict()["molecule_identifier"]
    new_name = inp.as_dict()["new_name"]

    if get_smol_from_mws(cmd_pointer, new_name, ignore_synonyms=True) is not None:
        output_error("There already exists a molecule named '{new_name}' in your working set.", return_val=False)
        return False

    openad_mol = get_smol_from_list(identifier, cmd_pointer.molecule_list)
    if openad_mol is not None:
        openad_mol["identifiers"]["name"] = new_name
        output_success("Molecule successfully renamed", return_val=False)
        return True

    output_error("No molecule '{identifier}' was found.", return_val=False)
    return False


def export_mws(cmd_pointer, inp):
    """
    Export molecule working set to a dataframe (Notebook) or CSV file (CLI).
    """

    if len(cmd_pointer.molecule_list) == 0:
        output_error("No molecules stored in your working set.", return_val=False)
        return True

    as_file = "file_name" in inp.as_dict()
    workspace_path = cmd_pointer.workspace_path()
    mws = cmd_pointer.molecule_list
    file_name = None
    ext = "csv"  # Default export file format

    # Return dataframe in Jupyter/API
    if GLOBAL_SETTINGS["display"] in ["notebook", "api"] and not as_file:
        return molset2dataframe(mws)

    # Save to file
    else:
        # File name provided
        if as_file:
            file_name = inp.as_dict()["file_name"]

        # No file name provided --> use default
        else:
            file_name = f"smol_export.{ext}"
            output_warning(msg("war_no_filename_provided", f"{file_name}.{ext}"), return_val=False, pad=0)

        # Detect file extension and strip it
        if file_name.lower().endswith(".molset.json"):
            file_name = file_name[:-12]
            ext = "molset.json"
        elif file_name.lower().endswith(".json"):
            file_name = file_name[:-5]
            ext = "molset.json"
        elif file_name.lower().endswith(".sdf"):
            file_name = file_name[:-4]
            ext = "sdf"
        elif file_name.lower().endswith(".csv"):
            file_name = file_name[:-4]
            ext = "csv"
        elif file_name.lower().endswith(".smi"):
            file_name = file_name[:-4]
            ext = "smi"

        # Find the next available file file path
        file_path = f"{workspace_path}/{file_name}.{ext}"
        if os.path.exists(file_path):
            i = 1
            while os.path.exists(f"{workspace_path}/{file_name}-{i}.{ext}"):
                i = i + 1
            file_path = f"{workspace_path}/{file_name}-{i}.{ext}"

        # Store the file
        if ext == "molset.json":
            success, err = save_molset_as_json(mws, file_path)
        elif ext == "sdf":
            success, err = save_molset_as_sdf(mws, file_path)
        elif ext == "csv":
            success, err = save_molset_as_csv(mws, file_path)
        elif ext == "smi":
            success, err = save_molset_as_smiles(mws, file_path)

        # Success
        if success:
            return output_success(
                f"Result set saved to workspace as <yellow>{file_path.split('/')[-1]}</yellow>", pad=0
            )

        # Error
        elif err:
            if "invalid_mols" in err:
                return output_error([err["error_msg"], err["invalid_mols"]])
            else:
                return output_error(err["error_msg"])


def property_retrieve(molecule_identifier, molecule_property, cmd_pointer):
    """retrieves a property for a molecule"""
    mol = get_smol_from_mws(cmd_pointer, molecule_identifier)
    if mol is None:
        mol = get_smol_from_pubchem(molecule_identifier)
        if mol is not None:
            return mol["properties"][molecule_property.lower()]

        else:
            output_error("molecule not available on pubchem", return_val=False)
            return None
    else:
        return mol["properties"][molecule_property.lower()]


def get_smol_prop(cmd_pointer, inp):
    """
    Get a property from a molecule.
    """

    result = None
    molecule_identifier = inp.as_dict()["molecule_identifier"]
    molecule_property = inp.as_dict()["property"]
    smol = find_smol(cmd_pointer, molecule_identifier)

    if not smol:
        return output_error("No molecule found for this identifier")

    # List of synonyms
    if molecule_property in PCY_IDFR.keys():
        result = smol["identifiers"][molecule_property]
    elif molecule_property == "synonyms":
        # if len(smol.get("synonyms", [])) > 0:
        if len(smol["synonyms"]) > 0:
            result = smol["synonyms"]
        else:
            return output_error("No synonyms found for this molecule")

    # Property

    elif molecule_property in smol["properties"]:
        result = smol["properties"][molecule_property]

    # Identifier
    elif molecule_property in smol["identifiers"]:
        result = smol["identifiers"][molecule_property]

    # Fail
    else:
        return output_error(
            f"We couldn't find the <yellow>{molecule_property}</yellow> property for the molecule <yellow>{molecule_identifier}</yellow>."
        )

    # Return
    if GLOBAL_SETTINGS["display"] == "api":
        return result
    if molecule_property == "synonyms":
        return list_columns(result)
    else:
        return output_text(result)


def get_smol_prop_lookup_error(cmd_pointer, inp):
    """
    Display an error message for an invalid molecule property.

    This is handled by a separate command function because
    the main command won't be reconized by pyparsing when
    the requested property is not in the property list.
    """
    output = "<error>The requested molecule property is not supported</error>"
    output = output + "\n\n<h1>Available properties:</h1>\n"
    output = output + list_columns(SMOL_PROPERTIES)
    return output_text(output)


# Launch molecule viewer and display molecule.
def show_mol(cmd_pointer, inp):
    from openad.gui.gui_launcher import gui_init

    molecule_identifier = inp.as_dict()["molecule_identifier"]
    path = "smol/" + urllib.parse.quote(molecule_identifier, safe="")
    gui_init(cmd_pointer, path)
    return True


# Launch molset viewer and display a molecule set file.
def show_molset(cmd_pointer, inp):
    from openad.gui.gui_launcher import gui_init

    molset_file = inp.as_dict()["molset_file"]

    path = "~/" + urllib.parse.quote(molset_file, safe="")
    gui_init(cmd_pointer, path)
    return True


# Launch molset viewer and display a molecule set dataframe.
def show_molset_df(cmd_pointer, inp):
    from openad.gui.gui_launcher import gui_init

    # molset_dataframe = cmd_pointer.api_variables[inp.as_dict()["in_dataframe"]]
    df_name = inp.as_dict()["in_dataframe"]

    path = "dataframe/" + df_name
    gui_init(cmd_pointer, path)
    return True


def _create_workspace_dir_if_nonexistent(cmd_pointer, dir_name):
    """creates a workspace directory"""
    if not os.path.isdir(cmd_pointer.workspace_path(cmd_pointer.settings["workspace"].upper()) + "/" + dir_name):
        os.mkdir(cmd_pointer.workspace_path(cmd_pointer.settings["workspace"].upper()) + "/" + dir_name)
    return cmd_pointer.workspace_path(cmd_pointer.settings["workspace"].upper()) + "/" + dir_name


#
#
# DEPRECATED
#
#


# Trash
# def _load_molecules(location):
#     """Loads molecules from  a given file"""
#     if not os.path.exists(os.path.expanduser(location)):
#         return None
#     with open(os.path.expanduser(location), "rb") as handle:
#         molecule = pickle.loads(handle.read())
#         handle.close()
#         return molecule
#


# MAJOR-RELEASE-TODO: Remove this, this is deprecated functionality
# We don't use .molecule files anymore
def save_molecules_DEPRECATED(cmd_pointer, inp):
    """saves a molecule set"""
    if "molset_name" not in inp:
        return False
    mol_file_path = _create_workspace_dir_if_nonexistent(cmd_pointer, "_mols")
    if cmd_pointer.molecule_list is not None and len(cmd_pointer.molecule_list) > 0:
        for mol in cmd_pointer.molecule_list:
            name = inp["molset_name"].upper() + "--" + mol["identifiers"]["inchikey"] + ".molecule"
            _write_molecules_DEPRECATED(mol, mol_file_path + "/" + name.strip())
    return True


# MAJOR-RELEASE-TODO: Remove this, this is deprecated functionality
# We don't use .molecule files anymore
def _write_molecules_DEPRECATED(molecule: dict, location):
    """writes molecules to a given file"""
    with open(os.path.expanduser(location), "wb") as handle:
        pickle.dump(molecule, handle)
    return True


# MAJOR-RELEASE-TODO: Remove this, this is deprecated functionality
# We don't use .molecule files anymore
def load_molecules_DEPRECATED(cmd_pointer, inp):
    """loads a molecule set into the molecule working set"""
    if "molset_name" not in inp:
        return False

    mol_file_path = _create_workspace_dir_if_nonexistent(cmd_pointer, "_mols")
    if "append" not in inp:
        cmd_pointer.molecule_list.clear()

    for i in glob.glob(mol_file_path + "/" + inp["molset_name"].upper() + "--*.molecule", recursive=True):
        func_file = open(i, "rb")
        mol = dict(pickle.load(func_file))
        for properties in mol["properties"]:
            if properties not in SMOL_PROPERTIES and properties not in [
                "atoms",
                "bonds",
                "record",
                "elements",
                "cactvs_fingerprint",
                "fingerprint",
            ]:
                SMOL_PROPERTIES.append(properties)
        cmd_pointer.molecule_list.append(mol.copy())
    output_text("<green>Number of molecules loaded</green> = " + str(len(cmd_pointer.molecule_list)), return_val=False)
    return True


# MAJOR-RELEASE-TODO: Remove this, this is deprecated functionality
# We don't use .molecule files anymore
def merge_molecules_DEPRECATED(cmd_pointer, inp):
    """loads a molecule set into the molecule working set"""
    if "molset_name" not in inp:
        return False
    merged = 0
    appended = 0
    mol_file_path = _create_workspace_dir_if_nonexistent(cmd_pointer, "_mols")

    for i in glob.glob(mol_file_path + "/" + inp["molset_name"].upper() + "--*.molecule", recursive=True):
        func_file = open(i, "rb")
        merge_mol = dict(pickle.load(func_file))
        existing_mol_twin = get_smol_from_mws(cmd_pointer, merge_mol["identifiers"]["canonical_smiles"])
        if existing_mol_twin is not None:
            if "append_only" not in inp.as_dict():
                merged += 1
                merge_molecule_REPLACE(merge_mol, existing_mol_twin)
        else:
            if "merge_only" not in inp.as_dict():
                cmd_pointer.molecule_list.append(merge_mol)
                appended += 1
    output_text("<green>Number of molecules added</green> = " + str(appended), return_val=False)
    output_text("<green>Number of molecules updated</green> = " + str(merged), return_val=False)
    return True


# MAJOR-RELEASE-TODO: Remove this, this is deprecated functionality
# We don't use .molecule files anymore
def display_molsets_DEPRECATED(cmd_pointer, inp):
    """displays the list of molecule-sets"""
    return _list_molsets_DEPRECATED(cmd_pointer)


# MAJOR-RELEASE-TODO: Remove this, this is deprecated functionality
# We don't use .molecule files anymore
def _list_molsets_DEPRECATED(cmd_pointer):
    """creates list of molecule sets"""
    molsets = []
    in_list = []
    mol_file_path = _create_workspace_dir_if_nonexistent(cmd_pointer, "_mols")
    for i in glob.glob(mol_file_path + "/*--*.molecule", recursive=True):
        x = i.split("/")
        molset = str(x[-1])

        molset = str(molset.split("--")[0])

        if molset not in in_list:
            in_list.append(molset)
            molsets.append([molset])
    if len(in_list) > 0:
        return output_table(molsets, is_data=False, headers=["Stored Molecule Sets"])
    return True
