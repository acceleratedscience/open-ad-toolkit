"""Functions that are called for molecule commands"""
import glob
import pickle
import sys, os
import shutil
import re
import json
import pandas as pd


from rdkit import Chem
from rdkit.Chem import AllChem


from openad.helpers.general import confirm_prompt
from openad.helpers.output import output_text, output_table, output_warning
from openad.molecules.mol_functions import canonical_smiles

from openad.molecules.mol_functions import (
    get_mol_from_formula,
    get_mol_from_inchi,
    get_mol_from_inchikey,
    get_mol_from_name,
    get_mol_from_smiles,
    get_mol_from_cid,
    new_molecule,
)
from openad.molecules.mol_functions import get_properties, get_identifiers

from openad.plugins.style_parser import print_s, style

cli_width = shutil.get_terminal_size().columns


class bold_style:
    """bold markers"""

    BOLD = "<b>"
    END = "</b>"


# from openad.molecules.rdkit_draw import print_mol_ascii


def display_molecule(cmd_pointer, inp):
    """displays a molecule and its properties"""
    molecule_identifier = inp.as_dict()["molecule_identifier"]

    mol = retrieve_mol_from_list(cmd_pointer, molecule_identifier)

    if mol is not None:
        print_string = format_identifers(mol) + "\n" + format_properties(mol) + "\n" + format_analysis(mol)

        # return print_s(print_string)
    else:
        mol = retrieve_mol(molecule_identifier)
        if mol is not None:
            cmd_pointer.last_external_molecule = mol.copy()
            print_string = format_identifers(mol) + "\n" + format_properties(mol) + "\n" + format_analysis(mol)
            # return print_s(print_string)
        else:
            print_s("molecule not available on pubchem")
            return None
    if cmd_pointer.notebook_mode is True:
        import py3Dmol
        from IPython.display import Markdown, display, HTML

        astyle = "stick"
        view_mol = Chem.MolFromSmiles(mol["properties"]["canonical_smiles"])  # pylint: disable=no-member
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
        print_string = print_string.replace("<success>", "<text style=color:green;white-space=pre>")
        print_string = print_string.replace("</success>", "</text>")
        print_string = print_string.replace("<yellow>", "<b>")
        print_string = print_string.replace("</yellow>", "</b>")
        print_string = print_string.replace("\n", "<br>")
        display(HTML("<pre>" + print_string + "</pre>"))
    else:
        print_s(print_string)

    return True


def export_molecule(cmd_pointer, inp):
    """exports a molecule as a dictionary"""
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
        print_s("file " + mol["name"] + ".json Saved to the current workspace")
    elif mol is not None and cmd_pointer.notebook_mode is True:
        return mol.copy()
    return True


def add_molecule(cmd_pointer, inp, force=False):
    """adds a molecule to the working set"""
    if isinstance(inp, dict):
        molecule_identifier = inp["molecule_identifier"]
    else:
        molecule_identifier = inp.as_dict()["molecule_identifier"]
    if (
        cmd_pointer.last_external_molecule is not None
        and is_molecule(cmd_pointer.last_external_molecule, molecule_identifier) is not None
    ):
        mol = cmd_pointer.last_external_molecule
    else:
        mol = retrieve_mol(molecule_identifier)
    if mol is None:
        output_text("Unable to identify molecule", cmd_pointer=cmd_pointer, return_val=False)
        return True

    identifier = mol["name"] + "   " + mol["properties"]["canonical_smiles"]

    if retrieve_mol_from_list(cmd_pointer, mol["properties"]["canonical_smiles"]) is not None:
        print_s("Molecule already in list")
        return True

    if force is False:
        if confirm_prompt("Are you wish to add " + identifier + " to your working list ?"):
            cmd_pointer.molecule_list.append(mol.copy())
            print_s("Molecule was Added.")
            return True

        print_s("Molecule was not added")
        return False
    cmd_pointer.molecule_list.append(mol.copy())
    return True


def create_molecule(cmd_pointer, inp):
    """creates a blank molecule"""

    molecule_identifier = inp.as_dict()["smiles"]
    molecule_name = inp.as_dict()["name"]

    mol = new_molecule(molecule_name, molecule_identifier)

    identifier = mol["name"] + "   " + mol["properties"]["canonical_smiles"]

    if retrieve_mol_from_list(cmd_pointer, mol["properties"]["canonical_smiles"]) != None:
        print_s("Molecule already in list")
        return True

    if confirm_prompt("Are you wish to add " + identifier + " to your working list ?"):
        cmd_pointer.molecule_list.append(mol.copy())
        print_s("Molecule was Added.")
        return True

    print_s("Molecule was not added")
    return False


def clear_workset(cmd_pointer, inp):
    """clears a workset"""
    if confirm_prompt("Are you wish to clear the Molecule Workset ?"):
        cmd_pointer.molecule_list.clear()


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
            print_s("Molecule was removed.")
        return True

    print_s("No Molecule Found")
    return True


def list_molecules(cmd_pointer, inp):
    """lists all molecules in the working set"""
    display_list = pd.DataFrame()

    if len(cmd_pointer.molecule_list) > 0:
        for mol in cmd_pointer.molecule_list:
            identifiers = get_identifiers(mol)
            display_list = pd.concat([display_list, pd.DataFrame([identifiers])])
        return display_list
    else:
        return output_text("No Molecules in List", cmd_pointer=cmd_pointer)


def retrieve_mol_from_list(cmd_pointer, molecule):
    """retrieves a molecule from the working list"""
    for mol in cmd_pointer.molecule_list:
        m = is_molecule(mol, molecule)
        if m is not None:
            return m.copy()

    for mol in cmd_pointer.molecule_list:
        m = is_molecule_synonym(mol, molecule)
        if m is not None:
            return m.copy()
    return None


def rename_mol_in_list(cmd_pointer, inp):
    """renames a molecule in the working list"""
    if retrieve_mol_from_list(cmd_pointer, inp.as_dict()["new_name"]) is not None:
        print_s("A molecule in Working Set already contains the new name.")
        return False

    for mol in cmd_pointer.molecule_list:
        m = is_molecule(mol, inp.as_dict()["molecule_identifier"])
        if m is not None:
            m["name"] = inp.as_dict()["new_name"]
            print_s("<success> molecule successfully re-named</success>")
            return True

    for mol in cmd_pointer.molecule_list:
        m = is_molecule_synonym(mol, inp.as_dict()["molecule_identifier"])
        if m is not None:
            m["name"] = inp.as_dict()["new_name"]
            print_s("<success> molecule successfully re-named</success>")
            return True

    print_s(" molecule was not renamed, no molecule found")

    return False


def export_molecule_set(cmd_pointer, inp):
    """exports molecule Set to Data frame on Notebook or file in workspace"""

    if len(cmd_pointer.molecule_list) == 0:
        print_s("\nNo Molecules in Molecule-Set")
        return True
    csv_file_name = None

    if cmd_pointer.notebook_mode and "csv_file_name" not in inp.as_dict():
        return moleculelist_to_data_frame(cmd_pointer.molecule_list.copy())
    else:
        if "csv_file_name" not in inp.as_dict():
            output_warning("WARNING no File Name Provided, reverting to Default.")
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
        print_s(f"Result set saved in Workspace as {file_name.split('/')[-1]}")


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


def is_molecule(mol, molecule):
    """determines if a given molecule identifier is actually a valid molecule"""
    if molecule.upper() == mol["name"].upper():
        return mol
    try:
        if int(molecule) == int(mol["properties"]["cid"]):
            return mol
    except:
        pass
    if molecule == mol["properties"]["inchi"]:
        return mol
    if molecule == mol["properties"]["inchikey"]:
        return mol
    if (
        mol["properties"]["isomeric_smiles"] != None
        and molecule.upper() == mol["properties"]["isomeric_smiles"].upper()
    ):
        return mol
    try:
        if canonical_smiles(molecule) == canonical_smiles(mol["properties"]["canonical_smiles"]):
            return mol
    except:
        pass
    return None


def is_molecule_synonym(mol, molecule):
    """determines if a molecule is mentioned in the synonym property of a molecule"""
    if mol["synonyms"] is not None and "Synonym" in mol["synonyms"]:
        for syn in mol["synonyms"]["Synonym"]:
            if molecule.upper() == syn.upper():
                return mol

    return None


def retrieve_mol(molecule):
    """gets molecule from pubchem"""
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
    """formats the identifiers for display"""
    id_string = "\n<yellow>Name:</yellow> {} \n".format(mol["name"])

    identifiers = get_identifiers(mol)
    i = 0
    for key, value in identifiers.items():
        if key != "name":
            if value is None:
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
    id_string = re.sub(r"<(.*?:)> ", r"<success>\1</success>", id_string)
    if mol["synonyms"] != None and "Synonym" in mol["synonyms"]:
        id_string = id_string + "\n\n<yellow>Synonyms:</yellow> {} \n\n".format(mol["synonyms"]["Synonym"])
    return id_string


def format_properties(mol):
    """formats properties for display"""
    id_string = "\n<yellow>Properties:</yellow>\n"
    identifiers = get_properties(mol)

    i = 0
    for key, value in identifiers.items():
        if key != "name":
            if value is None:
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
    id_string = re.sub(r"<(.*?:)> ", r"<success>\1</success>", id_string)
    return id_string


def format_analysis(mol):
    """formats analysis for display"""
    if "analysis" not in mol:
        return ""
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
    id_string = re.sub(r"<(.*?:)> ", r"<success>\1</success> ", id_string) + "\n"
    return id_string


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
    cmd_pointer.molecule_list.clear()

    for i in glob.glob(mol_file_path + "/" + inp["molecule-set_name"].upper() + "_*.molecule", recursive=True):
        func_file = open(i, "rb")
        mol = dict(pickle.load(func_file))
        cmd_pointer.molecule_list.append(mol.copy())
    print_s("\nNumber of Molecules Loaded = " + str(len(cmd_pointer.molecule_list)))
    return True


def display_molsets(cmd_pointer, inp):
    """displays the list of molecule-sets"""
    return list_molsets(cmd_pointer)


def list_molsets(cmd_pointer):
    """creates list of molecule sets"""
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
        return output_table(molsets, cmd_pointer=cmd_pointer, headers=["Stored Molecule Sets"])
    return True


def save_molecules(cmd_pointer, inp):
    """saves a molecule set"""
    if "molecule-set_name" not in inp:
        return False
    mol_file_path = _create_workspace_dir_if_nonexistent(cmd_pointer, "_mols")
    if cmd_pointer.molecule_list is not None and len(cmd_pointer.molecule_list) > 0:
        for mol in cmd_pointer.molecule_list:
            name = inp["molecule-set_name"].upper() + "_" + mol["properties"]["inchikey"] + ".molecule"
            _write_molecules(mol, mol_file_path + "/" + name.strip())
    return True


def property_retrieve(molecule_identifier, molecule_property, cmd_pointer):
    """retrieves a property for a molecule"""
    mol = retrieve_mol_from_list(cmd_pointer, molecule_identifier)
    if mol is None:
        mol = retrieve_mol(molecule_identifier)
        if mol is not None:
            return mol["properties"][molecule_property.lower()]

        else:
            print_s("molecule not available on pubchem")
            return None
    else:
        return mol["properties"][molecule_property.lower()]


def get_property(cmd_pointer, inp):
    """gets property from a molecule"""
    molecule_identifier = inp.as_dict()["molecule_identifier"]
    molecule_property = inp.as_dict()["property"]
    mol = retrieve_mol_from_list(cmd_pointer, molecule_identifier)
    if mol is None:
        mol = retrieve_mol(molecule_identifier)
        if mol is not None:
            # if not cmd_pointer.notebook_mode:
            #    print(mol["properties"][molecule_property.lower()])

            return mol["properties"][molecule_property.lower()]
            # return print_s(print_string)
        else:
            print_s("molecule not available on pubchem")
            return None
    else:
        return mol["properties"][molecule_property.lower()]


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
