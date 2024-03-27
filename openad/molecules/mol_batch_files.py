""" handles back load and unloading for molecule operations"""

import pandas
from rdkit.Chem import PandasTools
from openad.molecules.mol_functions import (
    merge_molecule_properties,
    valid_smiles,
    new_molecule,
    canonical_smiles,
    INPUT_MAPPINGS,
)
from openad.molecules.mol_commands import retrieve_mol_from_list, add_molecule
from openad.helpers.output import output_error, output_warning, output_success, output_text

from openad.helpers.output_msgs import msg
from openad.app.global_var_lib import GLOBAL_SETTINGS

naming_cache = {}


def load_batch_molecules(cmd_pointer, inp):
    """loads molecules in batch"""
    mol_dataframe = None
    if "load_molecules_dataframe" in inp.as_dict():
        mol_dataframe = _normalize_mol_df(cmd_pointer.api_variables[inp.as_dict()["in_dataframe"]], cmd_pointer)
    else:
        mol_dataframe = load_mol(inp.as_dict()["moles_file"], cmd_pointer)
    if mol_dataframe is None:
        output_error("Source not Found ", return_val=False)
        return True
    cmd_pointer.molecule_list.clear()
    if "pubchem_merge" in inp.as_dict():
        batch_pubchem(cmd_pointer, mol_dataframe)
    if mol_dataframe is not None:
        shred_merge_add_df_mols(mol_dataframe, cmd_pointer)
        output_success("Records loaded from file : " + str(len(mol_dataframe)), return_val=False)
        output_success(
            "Unique Molecules loaded to Working list : " + str(len(cmd_pointer.molecule_list)), return_val=False
        )

    return True


def batch_pubchem(cmd_pointer, dataframe):
    """does the prompting of pubchem for data to merge in a bach operation"""
    if GLOBAL_SETTINGS["display"] == "notebook":
        from halo import HaloNotebook as Halo  # pylint: disable=import-outside-toplevel
    else:
        from halo import Halo  # pylint: disable=import-outside-toplevel

    class Spinner(Halo):
        "contextual spinner"

        def __init__(self):
            # Alternative spinners:
            # simpleDotsScrolling, interval=100
            super().__init__(spinner="dots", color="white")

    batch_spinner = Spinner()
    batch_spinner.start("loading molecules from PubChem")
    dict_list = dataframe.to_dict("records")

    for a_mol in dict_list:
        try:
            Name_Flag = False
            if "name" in a_mol:
                name = a_mol["name"]
            elif "Name" in a_mol:
                name = a_mol["Name"]
            elif "NAME" in a_mol:
                name = a_mol["NAME"]
            elif "chemical_name" in a_mol:
                name = a_mol["chemical_name"]
            else:
                name = None
            if name is not None:
                merge_mol = retrieve_mol_from_list(cmd_pointer, name)
                if merge_mol is not None:
                    Name_Flag = True
            batch_spinner.text = "loading :" + a_mol["SMILES"]
            if not valid_smiles(a_mol["SMILES"]):
                output_warning(
                    "error merging SMILES: "
                    + a_mol["SMILES"]
                    + " Smiles has been excluded by load, it is not recognised by RDKIT",
                    return_val=False,
                )
                continue

            add_molecule(cmd_pointer, {"molecule_identifier": a_mol["SMILES"]}, force=True, suppress=True)

        except Exception as e:
            print(e)
            output_warning(
                "error merging SMILES: " + a_mol["SMILES"] + " Smiles has been excluded by load issues with input file",
                return_val=False,
            )

    batch_spinner.succeed("Finished loading from PubChem")
    batch_spinner.stop()


def shred_merge_add_df_mols(dataframe, cmd_pointer):
    """shreds the molecule relevent properties from dataframe and loads into molecules"""
    dict_list = dataframe.to_dict("records")
    merge_list = []
    for a_mol in dict_list:
        Name_Flag = False
        if "name" in a_mol:
            name = a_mol["name"]
        elif "Name" in a_mol:
            name = a_mol["Name"]
        elif "NAME" in a_mol:
            name = a_mol["NAME"]
        elif "chemical name" in a_mol:
            name = a_mol["chemical name"]
        elif "chemical_name" in a_mol:
            name = a_mol["chemical_name"]
        else:
            name = None
        if name == "":
            name = None
        if name is not None:
            merge_mol = retrieve_mol_from_list(cmd_pointer, name)
            if merge_mol is not None:
                Name_Flag = True
        if not valid_smiles(a_mol["SMILES"]):
            output_warning(
                "error merging SMILES: "
                + a_mol["SMILES"]
                + " Smiles has been excluded by load, it is not recognised by RDKIT",
                return_val=False,
            )
            continue
        merge_mol = retrieve_mol_from_list(cmd_pointer, a_mol["SMILES"])
        if Name_Flag is True and merge_mol is None:
            # output_error("There is already a molecule by the name, adding  increment to the name  " + name, return_val=False)
            i = 1
            while retrieve_mol_from_list(cmd_pointer, name + "-" + str(i)) is not None:
                i = i + 1
            name = name + "-" + str(i)

        # if Name_Flag is True and merge_mol["properties"]["canonical_smiles"] != canonical_smiles(a_mol["SMILES"]):
        #    output_error("There is already a molecule by the name, adding  increment to the name " + name, return_val=False)
        #    continue

        if merge_mol is None:
            if name is None:
                name = a_mol["SMILES"]
            merge_mol = new_molecule(a_mol["SMILES"], name)
        if merge_mol is not None:
            merge_mol = merge_molecule_properties(a_mol, merge_mol)
        else:
            output_warning(
                "error merging SMILES: "
                + a_mol["SMILES"]
                + " Smiles has been excluded by load, it is not recognised by RDKIT",
                return_val=False,
            )
            pass

        i = 0
        updated_flag = False

        while i < len(cmd_pointer.molecule_list):
            if (
                merge_mol["properties"]["canonical_smiles"]
                == cmd_pointer.molecule_list[i]["properties"]["canonical_smiles"]
            ):
                cmd_pointer.molecule_list[i] = merge_mol
                merge_list.append(merge_mol["name"])
                i = len(cmd_pointer.molecule_list)
                updated_flag = True
            i = i + 1
        if updated_flag is False:
            cmd_pointer.molecule_list.append(merge_mol)
    if len(merge_list) > 0:
        merge_list = list(set(merge_list))
        merge_list.sort()
        output_warning("The following molecules had 1 or more duplicate entries and were merged :", return_val=False)
        output_text("\n   - " + "\n   - ".join(merge_list), return_val=False)
    return True


def load_mol(source_file, cmd_pointer):
    """loads molecules from a souce file"""
    if source_file.split(".")[-1].lower() == "sdf":
        # From sdf file
        try:
            name = source_file.split("/")[-1]
            SDFFile = cmd_pointer.workspace_path(cmd_pointer.settings["workspace"].upper()) + "/" + name

            return _normalize_mol_df(PandasTools.LoadSDF(SDFFile), cmd_pointer)

        except BaseException as err:
            output_error(msg("err_load_sdf", err), return_val=False)
            return None
    elif source_file.split(".")[-1].lower() == "csv":
        # From csv file.
        try:
            name = source_file.split("/")[-1]
            SDFFile = cmd_pointer.workspace_path(cmd_pointer.settings["workspace"].upper()) + "/" + name

            mol_frame = pandas.read_csv(SDFFile)

            return _normalize_mol_df(mol_frame, cmd_pointer)

        except BaseException as err:
            output_error(msg("err_load_csv", err), return_val=False)
            return None


def _normalize_mol_df(mol_df: pandas.DataFrame, cmd_pointer):
    has_name = False
    contains_name = None

    for i in mol_df.columns:
        # Find the name column.
        if str(i.upper()) == "NAME" or str(i.lower()) == "chemical_name":
            has_name = True
        if contains_name is None and "NAME" in str(i.upper()):
            contains_name = i
        if contains_name is None and "CHEMICAL_NAME" in str(i.upper()):
            contains_name = i

        # Normalize any columns we'll be referring to later.
        if str(i.upper()) == "SMILES":
            mol_df.rename(columns={i: "SMILES"}, inplace=True)
        if str(i.upper()) == "ROMOL":
            mol_df.rename(columns={i: "ROMol"}, inplace=True)
        if str(i.upper()) == "IMG":
            mol_df.rename(columns={i: "IMG"}, inplace=True)
        if i in INPUT_MAPPINGS:
            mol_df.rename(columns={i: INPUT_MAPPINGS[i]}, inplace=True)

    # Normalize name column.
    if has_name is False and contains_name is not None:
        mol_df.rename(columns={contains_name: "NAME"}, inplace=True)

    # Add names when missing.
    try:
        if has_name is False:
            output_warning(msg("no_m2g_name_column"), return_val=False)

            mol_df["NAME"] = "unknown"
            for i in mol_df.itertuples():
                mol_df.loc[i.Index, "NAME"] = _smiles_to_iupac(mol_df.loc[i.Index, "SMILES"])

    except BaseException as err:
        return None
    return mol_df


def _smiles_to_iupac(smiles):
    import pubchempy

    if smiles in naming_cache:
        return naming_cache[smiles]
    try:
        compounds = pubchempy.get_compounds(smiles, namespace="smiles")
        match = compounds[0]
        naming_cache[smiles] = str(match)
    except BaseException:
        match = smiles
    return str(match)
