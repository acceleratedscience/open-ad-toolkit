""" handles back load and unloading for molecule operations"""

import pandas
from rdkit.Chem import PandasTools
from openad.molecules.mol_functions import (
    merge_molecule_properties,
    valid_smiles,
    new_molecule,
    canonical_smiles,
)
from openad.molecules.mol_commands import retrieve_mol_from_list, add_molecule
from openad.helpers.output import output_error, msg, output_warning
from openad.plugins.style_parser import print_s

naming_cache = {}


def load_batch_molecules(cmd_pointer, inp):
    """loads molecules in batch"""
    mol_dataframe = None
    if "load_molecules_dataframe" in inp.as_dict():
        mol_dataframe = cmd_pointer.api_variables[inp.as_dict()["in_dataframe"]]
    else:
        mol_dataframe = load_mol(inp.as_dict()["moles_file"], cmd_pointer)
    if mol_dataframe is None:
        print_s("\n Source not Found \n")
        return True
    if "pubchem_merge" in inp.as_dict():
        batch_pubchem(cmd_pointer, mol_dataframe)
    if mol_dataframe is not None:
        shred_merge_add_Dataframe_mols(mol_dataframe, cmd_pointer)
    return True


def batch_pubchem(cmd_pointer, dataframe):
    """does the prompting of pubchem for data to merge in a bach operation"""
    if cmd_pointer.notebook_mode is True:
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
            if "Name" in a_mol:
                name = a_mol["Name"]
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
                    cmd_pointer,
                    return_val=False,
                )
                continue

            add_molecule(cmd_pointer, {"molecule_identifier": a_mol["SMILES"]}, True)

        except Exception as e:
            print(e)
            output_warning(
                "error merging SMILES: " + a_mol["SMILES"] + " Smiles has been excluded by loadissues with input file",
                cmd_pointer,
                return_val=False,
            )

    batch_spinner.succeed("Finished loading from PubChem")
    batch_spinner.stop()


def shred_merge_add_Dataframe_mols(dataframe, cmd_pointer):
    """shreds the molecule relevent propoerties from data frame and loads into molecules"""
    dict_list = dataframe.to_dict("records")
    for a_mol in dict_list:
        Name_Flag = False
        if "name" in a_mol:
            name = a_mol["name"]
        if "Name" in a_mol:
            name = a_mol["Name"]
        else:
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
                cmd_pointer,
                return_val=False,
            )
            continue
        merge_mol = retrieve_mol_from_list(cmd_pointer, a_mol["SMILES"])
        if Name_Flag is True and merge_mol is None:
            print_s("There is already a molecule by the name " + name)
            continue
        if Name_Flag is True and merge_mol["properties"]["canonical_smiles"] != canonical_smiles(a_mol["SMILES"]):
            print_s("There is already a molecule by the name " + name)
            continue

        if merge_mol is None:
            if name is None:
                name = a_mol["SMILES"]
            merge_mol = new_molecule(name, a_mol["SMILES"])
        if merge_mol is not None:
            merge_mol = merge_molecule_properties(a_mol, merge_mol)
        else:
            output_warning(
                "error merging SMILES: "
                + a_mol["SMILES"]
                + " Smiles has been excluded by load, it is not recognised by RDKIT",
                cmd_pointer,
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
                i = len(cmd_pointer.molecule_list)
                updated_flag = True
            i = i + 1
        if updated_flag is False:
            cmd_pointer.molecule_list.append(merge_mol)

    return True


def load_mol(source_file, cmd_pointer):
    """loads molecules from a souce file"""
    if source_file.split(".")[-1].lower() == "sdf":
        # From sdf file
        try:
            name = source_file.split("/")[-1]
            SDFFile = cmd_pointer.workspace_path(cmd_pointer.settings["workspace"].upper()) + "/" + name

            return PandasTools.LoadSDF(SDFFile)

        except BaseException as err:
            output_error(msg("err_load_sdf", err, split=True), cmd_pointer, return_val=False)
            return None
    elif source_file.split(".")[-1].lower() == "csv":
        # From csv file.
        try:
            name = source_file.split("/")[-1]
            SDFFile = cmd_pointer.workspace_path(cmd_pointer.settings["workspace"].upper()) + "/" + name
            mol_frame = pandas.read_csv(SDFFile)
            return _normalize_mol_df(mol_frame, cmd_pointer)

        except BaseException as err:
            output_error(msg("err_load_csv", err, split=True), cmd_pointer, return_val=False)
            return None


def _normalize_mol_df(mol_df: pandas.DataFrame, cmd_pointer):
    has_name = False
    contains_name = None

    for i in mol_df.columns:
        # Find the name column.
        if str(i.upper()) == "NAME":
            has_name = True
        if contains_name is None and "NAME" in str(i.upper()):
            contains_name = i

        # Normalize any columns we'll be referring to later.
        if str(i.upper()) == "SMILES":
            mol_df.rename(columns={i: "SMILES"}, inplace=True)
        if str(i.upper()) == "ROMOL":
            mol_df.rename(columns={i: "ROMol"}, inplace=True)
        if str(i.upper()) == "IMG":
            mol_df.rename(columns={i: "IMG"}, inplace=True)

    # Normalize name column.
    if has_name is False and contains_name is not None:
        mol_df.rename(columns={contains_name: "NAME"}, inplace=True)

    # Add names when missing.
    try:
        if "NAME" not in mol_df.columns:
            output_warning(msg("no_m2g_name_column"), cmd_pointer)

            mol_df["NAME"] = "unknown"
            for i in mol_df.itertuples():
                # print(_smiles_to_iupac('CCC(COC(=O)CS)(C(=O)C(=O)CS)C(=O)C(=O)CS'))
                mol_df.loc[i.Index, "NAME"] = _smiles_to_iupac(mol_df.loc[i.Index, "SMILES"])

    except BaseException as err:
        print_s(err)
        return None


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
