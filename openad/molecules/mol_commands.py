"""Functions that are called for molecule commands"""
import glob
import pickle
import pandas as pd
import shutil
import re
import json

from openad.helpers.general import confirm_prompt
from openad.helpers.output import output_text, output_table
from openad.molecules.mol_functions import cannonical_smiles

cli_width = shutil.get_terminal_size().columns
from openad.molecules.mol_functions import (
    get_mol_from_formula,
    get_mol_from_inchi,
    get_mol_from_inchikey,
    get_mol_from_name,
    get_mol_from_smiles,
    get_mol_from_cid,
)
from openad.molecules.mol_functions import get_properties, get_identifiers
import sys, os
from openad.plugins.style_parser import print_s

# from openad.molecules.rdkit_draw import print_mol_ascii


def display_molecule(cmd_pointer, inp):
    molecule_identifier = inp.as_dict()["molecule_identifier"]
    mol = retrieve_mol_from_list(cmd_pointer, molecule_identifier)
    if mol is not None:
        print_string = format_identifers(mol) + "\n" + format_properties(mol) + "\n" + format_analysis(mol)
        return print_s(print_string)
    else:
        mol = retrieve_mol(molecule_identifier)
        if mol is not None:
            cmd_pointer.last_external_molecule = mol.copy()
            print_string = format_identifers(mol) + "\n" + format_properties(mol) + "\n" + format_analysis(mol)
            return print_s(print_string)
        else:
            print("molecule not available on pubchem")
        return None


def export_molecule(cmd_pointer, inp):
    molecule_identifier = inp.as_dict()["molecule_identifier"]
    mol = retrieve_mol_from_list(cmd_pointer, molecule_identifier)
    if mol is None:
        mol = retrieve_mol(molecule_identifier)
        if mol is not None:
            cmd_pointer.last_external_molecule = mol
    if mol is not None and ("as_file" in inp.as_dict() or cmd_pointer.notebook_mode is False):
        json_file = open(
            cmd_pointer.workspace_path(cmd_pointer.settings["workspace"].upper()) + "/" + mol["name"] + ".json",
            "w",
            encoding="utf-8",
        )
        json.dump(mol, json_file)
        print("file " + mol["name"] + ".json Saved to the current workspace")
    elif mol is not None and cmd_pointer.notebook_mode is True:
        return mol.copy()
    return True


def add_molecule(cmd_pointer, inp):
    molecule_identifier = inp.as_dict()["molecule_identifier"]
    if (
        cmd_pointer.last_external_molecule != None
        and is_molecule(cmd_pointer.last_external_molecule, molecule_identifier) != None
    ):
        mol = cmd_pointer.last_external_molecule
    else:
        mol = retrieve_mol(molecule_identifier)
    if mol == None:
        output_text("Unable to identify molecule", cmd_pointer=cmd_pointer, return_val=False)
        return True
    identifier = mol["name"] + "   " + mol["properties"]["canonical_smiles"]
    if retrieve_mol_from_list(cmd_pointer, mol["properties"]["canonical_smiles"]) != None:
        print("Molecule already in list")
        return True
    if confirm_prompt("Are you wish to add " + identifier + " to your working list ?"):
        cmd_pointer.molecule_list.append(mol.copy())
        print("Molecule was Added.")
        return True

    print("Molecule was not added")
    return False


def remove_molecule(cmd_pointer, inp):
    molecule_identifier = inp.as_dict()["molecule_identifier"]
    mol = retrieve_mol_from_list(cmd_pointer, molecule_identifier)

    if mol != None:
        identifier = mol["name"] + "   " + mol["properties"]["canonical_smiles"]

        if confirm_prompt("Are you wish to Remove " + identifier + " from your working list ?"):
            i = 0
            while (
                cmd_pointer.molecule_list[i]["properties"]["canonical_smiles"] != mol["properties"]["canonical_smiles"]
            ):
                i = i + 1

            cmd_pointer.molecule_list.pop(i)
            print("Molecule was removed.")
        return True

    print("No Molecule Found")
    return True


def list_molecules(cmd_pointer, inp):
    display_list = pd.DataFrame()

    if len(cmd_pointer.molecule_list) > 0:
        for mol in cmd_pointer.molecule_list:
            identifiers = get_identifiers(mol)
            display_list = pd.concat([display_list, pd.DataFrame([identifiers])])
        return display_list
    else:
        return output_text("No Molecules to List", cmd_pointer=cmd_pointer)


def retrieve_mol_from_list_old(cmd_pointer, molecule):
    for mol in cmd_pointer.molecule_list:
        if molecule.upper() == mol["name"].upper():
            return mol
        if molecule == mol["properties"]["cid"]:
            return mol
        if molecule == mol["properties"]["inchi"]:
            return mol
        if molecule == mol["properties"]["inchikey"]:
            return mol
        if molecule.upper() == mol["properties"]["isomeric_smiles"].upper():
            return mol
        if molecule.upper() == mol["properties"]["canonical_smiles"].upper():
            return mol
    return None


def retrieve_mol_from_list(cmd_pointer, molecule):
    for mol in cmd_pointer.molecule_list:
        m = is_molecule(mol, molecule)
        if m is not None:
            return m.copy()
    return None


def is_molecule(mol, molecule):
    if molecule.upper() == mol["name"].upper():
        return mol
    if molecule == mol["properties"]["cid"]:
        return mol
    if molecule == mol["properties"]["inchi"]:
        return mol
    if molecule == mol["properties"]["inchikey"]:
        return mol
    if molecule.upper() == mol["properties"]["isomeric_smiles"].upper():
        return mol
    try:
        if cannonical_smiles(molecule) == cannonical_smiles(mol["properties"]["canonical_smiles"]):
            return mol
    except:
        pass
    return None


def retrieve_mol(molecule):
    success, mol, comp = get_mol_from_name(molecule)
    if success:
        return mol
    success, mol, comp = get_mol_from_inchi(molecule)
    if success:
        return mol

    success, mol, comp = get_mol_from_smiles(molecule)
    if success:
        return mol
    # commented out until getting no time outs from pubchem
    # success, mol, comp = get_mol_from_formula(molecule)
    # if success:
    #    return mol

    success, mol, comp = get_mol_from_inchikey(molecule)
    if success:
        return mol

    success, mol, comp = get_mol_from_cid(molecule)
    if success:
        return mol
    return None


def format_identifers(mol):
    id_string = "\n<yellow>Name:</yellow> {} \n".format(mol["name"])
    identifiers = get_identifiers(mol)
    i = 0
    for key, value in identifiers.items():
        if key != "name":
            if len("{}: {} ".format(key, value)) < cli_width / 2 and i == 0:
                id_string = id_string + "{:<40}".format("<{}:> {} ".format(key, value))
                i = 1
            elif len(" {}: {} ".format(key, value)) < cli_width / 2 and i == 1:
                id_string = id_string + "{:<40}".format(" <{}:> {}".format(key, value)) + "\n"
                i = 0
            elif len("{}: {} ".format(key, value)) > cli_width / 2 and i == 1:
                id_string = id_string + "\n<{}:> {}".format(key, value) + "\n"
                i = 0
            elif len("{}: {} ".format(key, value)) > cli_width / 2 and i == 0:
                id_string = id_string + "{:<40}".format("<{}:> {}".format(key, value)) + "\n"
                i = 0
    id_string = re.sub(r"<(.*?:)> ", r"<green>\1</green> ", id_string)
    return id_string


def format_properties(mol):
    id_string = "\n<yellow>Properties:</yellow>\n"
    identifiers = get_properties(mol)

    i = 0
    for key, value in identifiers.items():
        if key != "name":
            if len("{}: {} ".format(key, value)) < cli_width / 2 and i == 0:
                id_string = id_string + "{:<40}".format("<{}:> {} ".format(key, value))
                i = 1
            elif len(" {}: {} ".format(key, value)) < cli_width / 2 and i == 1:
                id_string = id_string + "{:<40}".format(" <{}:> {}".format(key, value)) + "\n"
                i = 0
            elif len("{}: {} ".format(key, value)) > cli_width / 2 and i == 1:
                id_string = id_string + "\n<{}:> {}".format(key, value) + "\n"
                i = 0
            elif len("{}: {} ".format(key, value)) > cli_width / 2 and i == 0:
                id_string = id_string + "{:<40}".format("<{}:> {}".format(key, value)) + "\n"
                i = 0
    id_string = re.sub(r"<(.*?:)> ", r"<green>\1</green> ", id_string)
    return id_string


def format_analysis(mol):
    if mol["analysis"] == {}:
        return ""
    id_string = "\n<yellow>Analysis:</yellow>\n"
    i = 0

    for item in mol["analysis"]:
        id_string = (
            id_string
            + "\n<yellow>Toolkit: </yellow>"
            + item["toolkit"]
            + " <yellow>Function: </yellow>"
            + item["function"]
            + "\n"
        )
        for key, value in item.items():
            if key not in ["toolkit", "function"]:
                if value == None:
                    value = ""
                if len("{}: {} ".format(key, value)) < cli_width / 2 and i == 0:
                    id_string = id_string + "{:<40}".format("<{}:> {} ".format(key, value))
                    i = 1
                elif len(" {}: {} ".format(key, value)) < cli_width / 2 and i == 1:
                    id_string = id_string + "{:<40}".format(" <{}:> {}".format(key, value)) + "\n"
                    i = 0
                elif len("{}: {} ".format(key, value)) > cli_width / 2 and i == 1:
                    id_string = id_string + "\n<{}:> {}".format(key, value) + "\n"
                    i = 0
                elif len("{}: {} ".format(key, value)) > cli_width / 2 and i == 0:
                    id_string = id_string + "{:<40}".format("<{}:> {}".format(key, value)) + "\n"
                    i = 0
    id_string = re.sub(r"<(.*?:)> ", r"<green>\1</green> ", id_string) + "\n"
    return id_string


def _create_workspace_dir_if_nonexistent(cmd_pointer, dir_name):
    if not os.path.isdir(cmd_pointer.workspace_path(cmd_pointer.settings["workspace"].upper()) + "/" + dir_name):
        os.mkdir(cmd_pointer.workspace_path(cmd_pointer.settings["workspace"].upper()) + "/" + dir_name)
    return cmd_pointer.workspace_path(cmd_pointer.settings["workspace"].upper()) + "/" + dir_name


def load_molecules(cmd_pointer, inp):
    if "molecule-set_name" not in inp:
        return False
    mol_file_path = _create_workspace_dir_if_nonexistent(cmd_pointer, "_mols")
    cmd_pointer.molecule_list = []

    for i in glob.glob(mol_file_path + "/" + inp["molecule-set_name"].upper() + "_*.molecule", recursive=True):
        func_file = open(i, "rb")
        cmd_pointer.molecule_list.append(pickle.load(func_file))

    return True


def display_molsets(cmd_pointer, inp):
    return list_molsets(cmd_pointer)


def list_molsets(cmd_pointer):
    molsets = []
    in_list = []
    mol_file_path = _create_workspace_dir_if_nonexistent(cmd_pointer, "_mols")
    for i in glob.glob(mol_file_path + "/*.molecule", recursive=True):
        x = i.split("/")
        molset = str(x[-1])

        molset = str(molset.split("_")[0])

        if molset not in in_list:
            in_list.append(molset)
            molsets.append([molset])
    if len(in_list) > 0:
        output_table(molsets, cmd_pointer=cmd_pointer, headers=["Stored Molecule Sets"])
    return True


def save_molecules(cmd_pointer, inp):
    if "molecule-set_name" not in inp:
        return False
    mol_file_path = _create_workspace_dir_if_nonexistent(cmd_pointer, "_mols")
    if cmd_pointer.molecule_list is not None and len(cmd_pointer.molecule_list) > 0:
        for mol in cmd_pointer.molecule_list:
            name = inp["molecule-set_name"].upper() + "_" + mol["properties"]["inchikey"] + ".molecule"
            _write_molecules(mol, mol_file_path + "/" + name.strip())
    return True


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
