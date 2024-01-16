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

# Helpers
from openad.helpers.general import confirm_prompt
from openad.helpers.output import output_text, output_table, output_warning, output_error
from openad.helpers.output_msgs import msg
from openad.helpers.format_columns import single_value_columns, name_and_value_columns


# Molecule functions
from openad.molecules.mol_functions import (
    get_mol_from_formula,
    get_mol_from_inchi,
    get_mol_from_inchikey,
    get_mol_from_name,
    get_mol_from_smiles,
    get_mol_from_cid,
    new_molecule,
    get_properties,
    get_identifiers,
    canonical_smiles,
    get_mol_basic,
    mol2svg,
    mol2sdf,
)

# Globals
from openad.app.global_var_lib import GLOBAL_SETTINGS

# Flask
from openad.flask_apps import launcher
from openad.flask_apps.molviewer.routes import fetchRoutesMolViewer

# TEMP
from openad.plugins.style_parser import print_s, style

CLI_WIDTH = min(shutil.get_terminal_size().columns, 150)


class bold_style:
    """bold markers"""

    BOLD = "<b>"
    END = "</b>"


# from openad.molecules.rdkit_draw import print_mol_ascii


def display_molecule(cmd_pointer, inp):
    """displays a molecule and its properties"""
    molecule_identifier = inp.as_dict()["molecule_identifier"]
    if GLOBAL_SETTINGS["display"] == "notebook":
        global CLI_WIDTH
        CLI_WIDTH = 100

    mol = retrieve_mol_from_list(cmd_pointer, molecule_identifier)

    if mol is not None:
        print_string = (
            format_identifers(mol)
            + "\n"
            + format_synonyms(mol)
            + "\n"
            + format_properties(mol)
            + "\n"
            + format_analysis(mol)
        )

        # return print_s(print_string)
    else:
        mol = retrieve_mol(molecule_identifier)
        if mol is not None:
            cmd_pointer.last_external_molecule = mol.copy()
            print_string = (
                format_identifers(mol)
                + "\n"
                + format_synonyms(mol)
                + "\n"
                + format_properties(mol)
                + "\n"
                + format_analysis(mol)
            )
            # return print_s(print_string)
        else:
            output_error(msg("err_mol_not_on_pubchem"))
            return None
    if GLOBAL_SETTINGS["display"] == "notebook":
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
        output_text(print_string, edge=True)

    return True


def display_property_sources(cmd_pointer, inp):
    """displays a molecule properties sources"""
    if GLOBAL_SETTINGS["display"] == "notebook":
        global CLI_WIDTH
        CLI_WIDTH = 100
    molecule_identifier = inp.as_dict()["molecule_identifier"]

    mol = retrieve_mol_from_list(cmd_pointer, molecule_identifier)

    if mol is not None:
        print_string = format_sources(mol)

        # return print_s(print_string)
    else:
        mol = retrieve_mol(molecule_identifier)
        if mol is not None:
            cmd_pointer.last_external_molecule = mol.copy()
            print_string = format_sources(mol)
            # return print_s(print_string)
        else:
            print_s("molecule not available on pubchem")
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
    if mol is not None and ("as_file" in inp.as_dict() or GLOBAL_SETTINGS["display"] != "notebook"):
        json_file = open(
            cmd_pointer.workspace_path(cmd_pointer.settings["workspace"].upper()) + "/" + mol["name"] + ".json",
            "w",
            encoding="utf-8",
        )
        json.dump(mol, json_file)
        print_s("file " + mol["name"] + ".json Saved to the current workspace")
    elif mol is not None and GLOBAL_SETTINGS["display"] == "notebook":
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
        output_text("Unable to identify molecule", return_val=False)
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
        return output_table(display_list)
    else:
        return output_text("No Molecules in List")


def retrieve_mol_from_list(cmd_pointer, molecule):
    """retrieves a molecule from the working list"""

    for mol in cmd_pointer.molecule_list:
        m = is_molecule(mol, molecule)

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

    if GLOBAL_SETTINGS["display"] == "notebook" and "csv_file_name" not in inp.as_dict():
        return moleculelist_to_data_frame(cmd_pointer.molecule_list.copy())
    else:
        if "csv_file_name" not in inp.as_dict():
            output_warning(msg("war_no_filename_provided", "mols_export.csv"))
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
    id_string = id_string + name_and_value_columns(
        identifiers,
        cli_width=CLI_WIDTH,
        display_width=40,
        exclusions=["toolkit", "function"],
    )

    id_string = re.sub(r"<(.*?:)> ", r"<success>\1</success>", id_string)
    return id_string


def format_sources(mol):
    """formats the identifiers for display"""
    sources_string = "\n<yellow>Name:</yellow> {} \n".format(mol["name"])

    for mol_property, source in mol["property_sources"].items():
        if mol_property not in mol["properties"] or mol["properties"][mol_property] == None:
            continue
        sources_string = sources_string + "\n\n<yellow>Property:</yellow> {} \n".format(mol_property)
        sources_string = sources_string + name_and_value_columns(source, cli_width=CLI_WIDTH, display_width=30)

    sources = re.sub(r"<(.*?:)> ", r"<success>\1</success>", sources_string)
    return sources


def format_synonyms(mol):
    """formats synonyms for display"""
    synonyms_string = "\n<yellow>Synonyms:</yellow>\n"
    if "Synonym" in mol["synonyms"]:
        synonyms = mol["synonyms"]["Synonym"]
    else:
        synonyms = []
    all_synonyms = True
    new_synonyms = []
    for synonym in synonyms:
        if len(synonym) > 30:
            all_synonyms = False
            continue
        if len(str(synonym).strip()) == 0:
            continue
        new_synonyms.append(synonym)

    synonyms_string = synonyms_string + single_value_columns(new_synonyms, CLI_WIDTH, 30)
    name = mol["name"]
    if all_synonyms is False:
        synonyms_string = (
            synonyms_string
            + f"\n\n<soft>This list is truncated. To view all synonyms, run <cmd>@{name}>>synonyms</cmd></soft>\n"
        )
    synonyms_string = re.sub(r"<(.*?:)> ", r"<success>\1</success>", synonyms_string)
    return synonyms_string


def format_properties(mol):
    """formats properties for display"""
    properties_string = "\n<yellow>Properties:</yellow>\n"
    properties = get_properties(mol)

    properites_string = properties_string + name_and_value_columns(
        properties,
        cli_width=CLI_WIDTH,
        display_width=40,
        exclusions=["name", "canonical_smiles", "inchi", "inchikey", "cid", "isomeric_smiles"],
    )
    properties_string = re.sub(r"<(.*?:)> ", r"<success>\1</success>", properites_string)
    return properties_string


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
            + "\n\n<yellow>Toolkit: </yellow>"
            + item["toolkit"]
            + " <yellow>Function: </yellow>"
            + item["function"]
            + "\n"
        )
        id_string = id_string + name_and_value_columns(
            item, cli_width=CLI_WIDTH, display_width=50, exclusions=["toolkit", "function"], indent="    "
        )
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

    for i in glob.glob(mol_file_path + "/" + inp["molecule-set_name"].upper() + "--*.molecule", recursive=True):
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
            if molecule_property.lower() == "synonyms":
                if mol["synonyms"] is not None and "Synonym" in mol["synonyms"]:
                    return mol["synonyms"]["Synonym"]
            return mol["properties"][molecule_property.lower()]
        else:
            print_s("molecule not available on pubchem")
            return None
    else:
        if molecule_property.lower() == "synonyms":
            if mol["synonyms"] is not None and "Synonym" in mol["synonyms"]:
                return mol["synonyms"]["Synonym"]
        return mol["properties"][molecule_property.lower()]


# Launch molecule viewer.
def show_mol(cmd_pointer, inp):
    molecule_identifier = inp.as_dict()["molecule_identifier"]

    # Try loading the molecule from your working set.
    mol = retrieve_mol_from_list(cmd_pointer, molecule_identifier)

    # Try generating a basic molecule from RDKit.
    # Only works with InChI or SMILES as molecule_identifier.
    if mol is None:
        mol = get_mol_basic(molecule_identifier)

    # Fetch the molecule from PubChem,
    # The molecule_identifier is probably its name, CID or InChIKey.
    if mol is None:
        mol = retrieve_mol(molecule_identifier)

    # Render SVG and SDF
    mol_rdkit = Chem.MolFromInchi(mol["properties"]["inchi"])
    if mol_rdkit:
        mol_svg = mol2svg(mol_rdkit)
        mol_sdf = mol2sdf(mol_rdkit)
    else:
        mol_svg, mol_sdf = None, None

    # Load routes and launch browser UI.
    routes = fetchRoutesMolViewer(mol, mol_sdf, mol_svg)

    if GLOBAL_SETTINGS["display"] == "notebook":
        # Jupyter
        launcher.launch(cmd_pointer, routes, "molviewer")
    else:
        # CLI
        launcher.launch(cmd_pointer, routes, "molviewer")


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
