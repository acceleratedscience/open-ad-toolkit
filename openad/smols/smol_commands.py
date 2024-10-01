"""Functions that are called for molecule commands"""

import glob
import pickle
import os
import json
import urllib.parse
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

# Helpers
from openad.helpers.general import confirm_prompt
from openad.helpers.output import output_text, output_table, output_warning, output_error, output_success
from openad.helpers.output_msgs import msg
from openad.helpers.pretty_data import list_columns, key_val_columns, key_val_full
from openad.app.global_var_lib import GLOBAL_SETTINGS


# Molecule functions
from openad.smols.smol_functions import (
    SMOL_PROPERTIES,
    find_smol,
    get_smol_from_mws,
    get_smol_from_pubchem,
    get_human_properties,
    _get_identifiers,
    get_smol_from_list,
    mws_add,
    mws_remove,
    merge_molecule_REPLACE,
)


def display_molecule(cmd_pointer, inp):
    """
    Display a molecule's identifiers, synonyms and properties.
    """

    molecule_identifier = inp.as_dict()["molecule_identifier"]

    # if GLOBAL_SETTINGS["display"] == "notebook":
    #     global PRINT_WIDTH
    #     PRINT_WIDTH = 100

    smol = find_smol(cmd_pointer, molecule_identifier, show_spinner=True)
    if smol is None:
        output_error("Error: Not a valid Molecule", return_val=False)
        return None

    if smol is not None:
        ouput = format_identifiers(smol) + "\n\n" + format_synonyms(smol) + "\n\n" + format_properties(smol)
        if "analysis" in smol and smol["analysis"] != []:
            ouput = ouput + "\n\n" + format_analysis(smol)

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
        ouput = ouput.replace("<soft>", "<span style=color:#cccccc;white-space=pre>")
        ouput = ouput.replace("</soft>", "</span>")
        ouput = ouput.replace("<cyan>", "<span style=color:#00AAAA;white-space=pre>")
        ouput = ouput.replace("</cyan>", "</span>")
        ouput = ouput.replace("\n", "<br>")
        display(HTML("<pre>" + ouput + "</pre>"))

        # display(output_text(print_string, return_val=True))
    else:
        output_text(ouput, edge=True, pad=1)

    return True


def format_identifiers(smol):
    """
    Format the identifiers for display.
    """

    output = f"<h1>{smol['name'].capitalize()}</h1>"
    identifiers = smol.get("identifiers", {})
    output = output + key_val_full(identifiers, print_width=GLOBAL_SETTINGS["print_width"] - 5)
    return output


def format_synonyms(smol, truncate=30):
    """
    Format synonyms for display.
    """

    synonyms = smol.get("synonyms", [])
    synonym_count = len(synonyms)

    # Truncate list of synonyms
    if truncate is not None:
        synonyms = synonyms[:truncate]
        is_truncated = synonym_count > len(synonyms)
    else:
        is_truncated = False

    # Title
    output = "<h1>Synonyms:</h1>"

    # Truncation note
    if is_truncated:
        output = (
            output
            + f"\n<soft>This list is truncated. To view all synonyms, run <cmd>@{smol['name']}>>synonyms</cmd></soft>\n\n"
        )
    else:
        output = output + "\n"

    # Values
    output = output + list_columns(synonyms, print_width=GLOBAL_SETTINGS["print_width"] - 5, is_truncated=is_truncated)

    return output


def format_properties(smol):
    """
    Format properties for display.
    """

    output = "<h1>Properties:</h1>"
    properties = get_human_properties(smol)
    output = output + key_val_columns(
        properties, print_width=GLOBAL_SETTINGS["print_width"] - 5, ignore_keys=["DS_URL"]
    )
    return output


def format_analysis(smol):
    """
    Format analysis for display.
    """

    if "analysis" not in smol or smol["analysis"] == []:
        return ""

    output = "<h1>Analysis:</h1>"

    for item in smol["analysis"]:
        output = (
            output
            + "\n<cyan>Toolkit: </cyan>"
            + item["toolkit"]
            + " / <cyan>Function: </cyan>"
            + item["function"]
            + "\n"
        )
        output = output + key_val_columns(item, print_width=GLOBAL_SETTINGS["print_width"] - 5, indent=4)
    return output


def display_property_sources(cmd_pointer, inp):
    """displays a molecule properties sources"""
    # if GLOBAL_SETTINGS["display"] == "notebook":
    #     global PRINT_WIDTH
    #     PRINT_WIDTH = 100
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


def format_sources(mol):
    """
    Format the sources for display.
    """

    output = f"\n<yellow>Name:</yellow> {mol['name']}\n"

    for mol_property, source in mol["property_sources"].items():
        if mol_property not in mol["properties"] or mol["properties"][mol_property] is None:
            continue
        output = output + "\n\n<yellow>Property:</yellow> {} \n".format(mol_property)
        output = output + key_val_columns(source, print_width=GLOBAL_SETTINGS["print_width"] - 5)

    return output


def export_molecule(cmd_pointer, inp):
    """exports a molecule as a dictionary"""

    molecule_identifier = inp.as_dict()["molecule_identifier"]
    mol = get_smol_from_mws(cmd_pointer, molecule_identifier)
    if mol is None:
        mol = get_smol_from_pubchem(molecule_identifier)
        if mol is not None:
            cmd_pointer.last_external_molecule = mol
    if mol is not None and ("as_file" in inp.as_dict() or GLOBAL_SETTINGS["display"] not in ["api", "notebook"]):
        json_file = open(
            cmd_pointer.workspace_path(cmd_pointer.settings["workspace"].upper()) + "/" + mol["name"] + ".json",
            "w",
            encoding="utf-8",
        )
        json.dump(mol, json_file)
        output_success("File " + mol["name"] + ".json saved to the current workspace", return_val=False)
    elif mol is not None and GLOBAL_SETTINGS["display"] in ["api", "notebook"]:
        return mol.copy()
    return True


def add_molecule(cmd_pointer, inp):
    """
    Adds a molecule to your molecules working set (my-mols).
    """

    identifier = inp.as_dict()["molecule_identifier"]

    if "basic" in inp.as_dict():
        basic = True
    else:
        basic = False

    if "name" in inp.as_dict():
        name = inp.as_dict()["name"]
    else:
        name = identifier

    if "force" in inp.as_dict():
        force = True
    else:
        force = False

    # Create molecule dict.
    openad_mol = find_smol(cmd_pointer, identifier, name, basic)
    # Add it to the working set.
    mws_add(cmd_pointer, openad_mol, force=force)


def remove_molecule(cmd_pointer, inp):
    """
    Removes a molecule from your molecules working set (my-mols).
    """

    if "force" in inp.as_dict():
        force = True
    else:
        force = False

    molecule_identifier = inp.as_dict()["molecule_identifier"]
    mol = get_smol_from_mws(cmd_pointer, molecule_identifier)
    mws_remove(cmd_pointer, mol, force=force)


# def create_molecule(cmd_pointer, inp):
#     """creates a blank molecule"""

#     molecule_identifier = inp.as_dict()["smiles"]
#     molecule_name = inp.as_dict()["name"]

#     mol = new_molecule(molecule_name, name=molecule_identifier)

#     identifier = mol["name"] + "   " + mol["properties"]["canonical_smiles"]

#     if get_smol_from_mws(cmd_pointer, mol["properties"]["canonical_smiles"]) != None:
#         output_error("Molecule already in list", return_val=False)
#         return True

#     if confirm_prompt("Are you wish to add " + identifier + " to your working list ?"):
#         cmd_pointer.molecule_list.append(mol.copy())
#         output_success("Molecule was added.", return_val=False)
#         return True

#     output_error("Molecule was not added", return_val=False)
#     return False


def clear_workset(cmd_pointer, inp):
    """clears a workset"""
    if confirm_prompt("Are you wish to clear the Molecule Workset ?"):
        cmd_pointer.molecule_list.clear()


def list_molecules(cmd_pointer, inp):
    """lists all molecules in the working set"""
    display_list = pd.DataFrame()

    if len(cmd_pointer.molecule_list) > 0:
        for mol in cmd_pointer.molecule_list:
            identifiers = _get_identifiers(mol)
            display_list = pd.concat([display_list, pd.DataFrame([identifiers])])
        return display_list
        # if GLOBAL_SETTINGS["display"] == "notebook":
        #    return output_table(display_list, is_data=True)
        # else:
        #    return output_table(display_list)

    else:
        return output_warning("Your molecules working set is empty")


def show_molecules(cmd_pointer, inp):
    """Display your working set in the browser or iframe"""
    from openad.gui.gui_launcher import gui_init

    gui_init(cmd_pointer, "my-mols")


def rename_mol_in_list(cmd_pointer, inp):
    """
    Renames a molecule in your working molecule list.
    """
    if get_smol_from_mws(cmd_pointer, inp.as_dict()["new_name"], ignore_synonyms=True) is not None:
        output_error("A molecule in your working set already contains the new name", return_val=False)
        return False

    identifier = inp.as_dict()["molecule_identifier"]
    openad_mol = get_smol_from_list(identifier, cmd_pointer.molecule_list)
    if openad_mol is not None:
        openad_mol["name"] = inp.as_dict()["new_name"]
        output_success("Molecule successfully renamed", return_val=False)
        return True

    # TRASH
    # @refactored
    # for mol in cmd_pointer.molecule_list:
    #     m = is_molecule(mol, inp.as_dict()["molecule_identifier"])
    #     if m is not None:
    #         m["name"] = inp.as_dict()["new_name"]
    #         output_success("Molecule successfully renamed", return_val=False)
    #         return True

    # TRASH
    # This is not integrated into mol_functions --> get_smol_from_list
    # for mol in cmd_pointer.molecule_list:
    #     m = is_molecule_synonym(mol, inp.as_dict()["molecule_identifier"])
    #     if m is not None:
    #         m["name"] = inp.as_dict()["new_name"]
    #         output_success("molecule successfully renamed", return_val=False)
    #         return True

    output_error("Molecule was not renamed, no molecule found", return_val=False)

    return False


def export_molecule_set(cmd_pointer, inp):
    """exports molecule Set to Data frame on Notebook or file in workspace"""

    if len(cmd_pointer.molecule_list) == 0:
        output_error("No molecules in molecule-set", return_val=False)
        return True
    csv_file_name = None

    if GLOBAL_SETTINGS["display"] in ["notebook", "api"] and "csv_file_name" not in inp.as_dict():
        return moleculelist_to_data_frame(cmd_pointer.molecule_list.copy())

    else:
        if "csv_file_name" not in inp.as_dict():
            output_warning(msg("war_no_filename_provided", "mols_export.csv"), return_val=False)
            csv_file_name = "mols_export"
        else:
            csv_file_name = inp.as_dict()["csv_file_name"]
        if csv_file_name.lower().endswith(".csv"):
            csv_file_name = csv_file_name.split(".")[0]

        file_name = cmd_pointer.workspace_path(cmd_pointer.settings["workspace"].upper()) + "/" + csv_file_name + ".csv"

        if os.path.exists(file_name):
            i = 0
            while os.path.exists(
                cmd_pointer.workspace_path(cmd_pointer.settings["workspace"].upper())
                + "/"
                + csv_file_name
                + str(i)
                + ".csv"
            ):
                i = i + 1
            file_name = (
                cmd_pointer.workspace_path(cmd_pointer.settings["workspace"].upper())
                + "/"
                + csv_file_name
                + str(i)
                + ".csv"
            )
        result = moleculelist_to_data_frame(cmd_pointer.molecule_list.copy())
        result.to_csv(file_name)
        output_success(f"Result set saved to workspace as {file_name.split('/')[-1]}", return_val=False)


def moleculelist_to_data_frame(molecule_set):
    """Turns the molecule list properties to a dataframe"""
    molecule_list = []
    for molecule in molecule_set:
        mol = {"mol_name": molecule["name"]}
        mol["SMILES"] = molecule["properties"]["canonical_smiles"]
        mol["canonical_smiles"] = molecule["properties"]["canonical_smiles"]
        mol["isomeric_smiles"] = molecule["properties"]["isomeric_smiles"]
        mol["inchi"] = molecule["properties"]["inchi"]
        mol["inchikey"] = molecule["properties"]["inchikey"]
        mol["formula"] = molecule["properties"]["molecular_formula"]

        for mol_property in molecule["properties"]:
            if mol_property not in mol:
                mol[mol_property] = molecule["properties"][mol_property]
        molecule_list.append(mol.copy())
    return pd.DataFrame(molecule_list)


# TRASH
# @refactored
# Moved to mol_functions --> get_smol_from_list (includes the loop)
# def is_molecule(mol, molecule):
#     """determines if a given molecule identifier is actually a valid molecule"""
#     if molecule.upper() == mol["name"].upper():
#         return mol
#     try:
#         if int(molecule) == int(mol["properties"]["cid"]):
#             return mol
#     except:
#         pass
#     if molecule == mol["properties"]["inchi"]:
#         return mol
#     if molecule == mol["properties"]["inchikey"]:
#         return mol
#     if (
#         mol["properties"]["isomeric_smiles"] is not None
#         and molecule.upper() == mol["properties"]["isomeric_smiles"].upper()
#     ):
#         return mol
#     try:
#         if canonical_smiles(molecule) == canonical_smiles(mol["properties"]["canonical_smiles"]):
#             return mol
#     except:
#         pass
#     return None


# TRASH
# @refactored
# This is not integrated into mol_functions --> get_smol_from_list
# def is_molecule_synonym(mol, molecule):
#     """determines if a molecule is mentioned in the synonym property of a molecule"""
#     if mol["synonyms"] is not None and "Synonym" in mol["synonyms"]:
#         for syn in mol["synonyms"]["Synonym"]:
#             if molecule.upper() == syn.upper():
#                 return mol

#     return None


def _create_workspace_dir_if_nonexistent(cmd_pointer, dir_name):
    """creates a workspace directory"""
    if not os.path.isdir(cmd_pointer.workspace_path(cmd_pointer.settings["workspace"].upper()) + "/" + dir_name):
        os.mkdir(cmd_pointer.workspace_path(cmd_pointer.settings["workspace"].upper()) + "/" + dir_name)
    return cmd_pointer.workspace_path(cmd_pointer.settings["workspace"].upper()) + "/" + dir_name


def load_molecules(cmd_pointer, inp):
    """loads a molecule set into the working list"""
    if "molecule-set_name" not in inp:
        return False

    mol_file_path = _create_workspace_dir_if_nonexistent(cmd_pointer, "_mols")
    if "append" not in inp:
        cmd_pointer.molecule_list.clear()

    for i in glob.glob(mol_file_path + "/" + inp["molecule-set_name"].upper() + "--*.molecule", recursive=True):
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


def display_molsets(cmd_pointer, inp):
    """displays the list of molecule-sets"""
    return list_molsets(cmd_pointer)


def list_molsets(cmd_pointer):
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


def save_molecules(cmd_pointer, inp):
    """saves a molecule set"""
    if "molecule-set_name" not in inp:
        return False
    mol_file_path = _create_workspace_dir_if_nonexistent(cmd_pointer, "_mols")
    if cmd_pointer.molecule_list is not None and len(cmd_pointer.molecule_list) > 0:
        for mol in cmd_pointer.molecule_list:
            name = inp["molecule-set_name"].upper() + "--" + mol["properties"]["inchikey"] + ".molecule"
            _write_molecules(mol, mol_file_path + "/" + name.strip())
    return True


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
    molecule_property = inp.as_dict()["property"].lower()
    smol = find_smol(cmd_pointer, molecule_identifier)

    if not smol:
        return output_error("No molecule found for this identifier")

    # List of synonyms
    if molecule_property == "synonyms":
        if len(smol.get("synonyms", [])) > 0:
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
        output_error(
            f"We couldn't find the '<reset>{molecule_property}</reset>' property for the molecule <yellow>{molecule_identifier}</yellow>."
        )

    # Return
    if GLOBAL_SETTINGS["display"] == "api":
        return result
    if molecule_property == "synonyms":
        return output_text(list_columns(result), pad=1)
    else:
        return output_text(result)


# Launch molecule viewer and display molecule.
def show_mol(cmd_pointer, inp):
    from openad.gui.gui_launcher import gui_init

    molecule_identifier = inp.as_dict()["molecule_identifier"]
    path = "smol/" + urllib.parse.quote(molecule_identifier, safe="")
    gui_init(cmd_pointer, path)


# Launch molset viewer and display a molecule set file.
def show_molset(cmd_pointer, inp):
    from openad.gui.gui_launcher import gui_init

    molset_file = inp.as_dict()["molset_file"]

    path = "~/" + urllib.parse.quote(molset_file, safe="")
    gui_init(cmd_pointer, path)


# Launch molset viewer and display a molecule set dataframe.
def show_molset_df(cmd_pointer, inp):
    from openad.gui.gui_launcher import gui_init

    # molset_dataframe = cmd_pointer.api_variables[inp.as_dict()["in_dataframe"]]
    df_name = inp.as_dict()["in_dataframe"]

    path = "dataframe/" + df_name
    gui_init(cmd_pointer, path)


def _load_molecules(location):
    """Loads molecules from  a given file"""
    if not os.path.exists(os.path.expanduser(location)):
        return None
    with open(os.path.expanduser(location), "rb") as handle:
        molecule = pickle.loads(handle.read())
        handle.close()
        return molecule


def _write_molecules(molecule: dict, location):
    """writes molecules to a given file"""
    with open(os.path.expanduser(location), "wb") as handle:
        pickle.dump(molecule, handle)
    return True


def merge_molecules(cmd_pointer, inp):
    """loads a molecule set into the working list"""
    if "molecule-set_name" not in inp:
        return False
    merged = 0
    appended = 0
    mol_file_path = _create_workspace_dir_if_nonexistent(cmd_pointer, "_mols")

    for i in glob.glob(mol_file_path + "/" + inp["molecule-set_name"].upper() + "--*.molecule", recursive=True):
        func_file = open(i, "rb")
        merge_mol = dict(pickle.load(func_file))
        existing_mol_twin = get_smol_from_mws(cmd_pointer, merge_mol["properties"]["canonical_smiles"])
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
