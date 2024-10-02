""" handles back load and unloading for molecule operations"""

import pandas
from rdkit import RDLogger
from rdkit.Chem import PandasTools
from openad.smols.smol_functions import (
    find_smol,
    get_smol_from_mws,
    merge_molecule_properties,
    valid_smiles,
    new_smol,
    mws_add,
    normalize_mol_df,
    canonicalize,
    load_mols_from_file,
)
from openad.smols.smol_transformers import dataframe2molset
from openad.app.global_var_lib import GLOBAL_SETTINGS
from openad.helpers.output import output_error, output_warning, output_success, output_text
from openad.helpers.output_msgs import msg
from openad.helpers.spinner import Spinner

RDLogger.DisableLog("rdApp.error")  # Suppress RDKiot errors

mol_name_cache = {}


#
#


def merge_molecule_property_data(cmd_pointer, inp):
    "merges data where SMILES,Property and Value are in Data Frame or csv"
    mol_dataframe = None
    if "force" in inp.as_dict():
        force = True
    else:
        force = False

    if "merge_molecules_data_dataframe" in inp.as_dict():
        mol_dataframe = cmd_pointer.api_variables[inp.as_dict()["in_dataframe"]]
    else:
        mol_dataframe = load_mol_data(inp.as_dict()["moles_file"], cmd_pointer)
    if mol_dataframe is None:
        output_error("Source not Found ", return_val=False)
        return True

    if "subject" in mol_dataframe.columns:
        SMILES = "subject"
    elif "smiles" in mol_dataframe.columns:
        SMILES = "smiles"
    elif "SMILES" in mol_dataframe.columns:
        SMILES = "SMILES"
    else:
        output_error("No  'subject' or 'SMILES' column found ", return_val=False)
        return True
    if "property" in mol_dataframe.columns:
        prop = "property"
    elif "PROPERTY" in mol_dataframe.columns:
        prop = "PROPERTY"
    else:
        output_error("No  'property' or 'PROPERTY' column found ", return_val=False)
        return True
    if "value" in mol_dataframe.columns:
        val = "value"
    elif "VALUE" in mol_dataframe.columns:
        val = "VALUE"
    elif "result" in mol_dataframe.columns:
        val = "result"
    elif "RESULT" in mol_dataframe.columns:
        val = "RESULT"
    else:
        output_error("No  'result' or 'value' column found ", return_val=False)
        return True

    mol_dataframe = mol_dataframe.pivot_table(index=SMILES, columns=[prop], values=val, aggfunc="first")

    mol_dataframe = mol_dataframe.reset_index()
    for row in mol_dataframe.to_dict("records"):
        update_flag = True
        merge_mol = None
        try:
            smiles = canonicalize(row[SMILES])

            merge_mol = get_smol_from_mws(cmd_pointer, smiles)
        except:
            output_warning("unable to canonicalise:" + row[SMILES])
            continue
        if merge_mol is None:
            merge_mol = new_smol(smiles, name=row[SMILES])

            update_flag = False
        else:
            update_flag = True

        # else duplicate
        # a_mol = {"SMILES": row[SMILES], row[prop]: row[val]}
        if merge_mol is not None:
            merge_mol = merge_molecule_properties(row, merge_mol)
            # print("updated: " + str(a_mol))
            if update_flag is False:
                cmd_pointer.molecule_list.append(merge_mol)
    output_success("Molecule Data Merged", return_val=False)
    return True


def load_mol_data(source_file, cmd_pointer):
    """loads molecule data from a file where Smiles, property and values are supplied in row format"""

    # SDF
    if source_file.split(".")[-1].lower() == "sdf":
        try:
            name = source_file.split("/")[-1]
            sdffile = cmd_pointer.workspace_path(cmd_pointer.settings["workspace"].upper()) + "/" + name
            mol_frame = PandasTools.LoadSDF(sdffile)
        except Exception as err:
            output_error(msg("err_load", "SDF", err), return_val=False)
            return None

    # CSV
    elif source_file.split(".")[-1].lower() == "csv":
        try:
            name = source_file.split("/")[-1]
            name = cmd_pointer.workspace_path(cmd_pointer.settings["workspace"].upper()) + "/" + name
            mol_frame = pandas.read_csv(name, dtype="string")
        except Exception as err:
            output_error(msg("err_load", "CSV", err), return_val=False)
            return None

    return mol_frame


def load_mols_to_mws(cmd_pointer, inp):
    """
    Load a batch of molecules into the molecules working set.
    """

    molset = None
    df_name = inp.as_dict().get("in_dataframe", None)
    file_path = inp.as_dict().get("moles_file", None)

    # Load from dataframe
    if df_name:
        df = cmd_pointer.api_variables[df_name]
        molset = dataframe2molset(df)
        # molset = normalize_mol_df(cmd_pointer.api_variables[inp.as_dict()["in_dataframe"]], batch=True)
        if molset is None:
            return output_error("The provided dataframe does not contain molecules")

    # Load from file
    elif file_path:
        molset = load_mols_from_file(cmd_pointer, file_path)
        if molset is None:
            return output_error("Source not Found")

    # Add molecules
    if "append" in inp:
        cmd_pointer.molecule_list.extend(molset)
    else:
        cmd_pointer.molecule_list = molset

    if "pubchem_merge" in inp.as_dict():
        batch_pubchem(cmd_pointer, molset)

    if molset:
        # shred_merge_add_df_mols(molset, cmd_pointer)
        # Todo: `load mols using file` should add instead of overwrite your current mols,
        # when this is implemented, we'll need to calculate successfully loaded mols differently.
        return output_success(
            f"Successfully loaded <yellow>{len(molset)}</yellow> molecules into the working set", pad=0
        )


def batch_pubchem(cmd_pointer, dataframe):
    """does the prompting of pubchem for data to merge in a bach operation"""

    if GLOBAL_SETTINGS["display"] == "notebook":
        from halo import HaloNotebook as Halo  # pylint: disable=import-outside-toplevel
    else:
        from halo import Halo  # pylint: disable=import-outside-toplevel

    batch_spinner = Spinner(GLOBAL_SETTINGS["verbose"])

    print(999, dataframe)

    batch_spinner.start("loading molecules from PubChem")

    dict_list = dataframe.to_dict("records")

    for i, a_mol in enumerate(dict_list):
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
                merge_mol = get_smol_from_mws(cmd_pointer, name)
                if merge_mol is not None:
                    Name_Flag = True
            batch_spinner.text = f"Loading: {a_mol['SMILES']}"

            if not valid_smiles(a_mol["SMILES"]):
                output_warning(
                    "error merging SMILES: "
                    + a_mol["SMILES"]
                    + " Smiles has been excluded by load, it is not recognised by RDKIT",
                    return_val=False,
                )
                continue

            # TRASH - @refactored
            # add_molecule(cmd_pointer, {"molecule_identifier": a_mol["SMILES"]}, force=True, suppress=True)

            # Create molecule dict.
            openad_mol = find_smol(cmd_pointer, a_mol["SMILES"])  # @@TODO: to be tested

            # Add it to the working set.
            mws_add(cmd_pointer, openad_mol, force=True, suppress=True)

        except Exception as err:
            print(err)
            err_msg = f"#{i} - <error>Invalid SMILES, molecule discarded:</error> <yellow>{a_mol['SMILES']}</yellow>"
            output_text(err_msg, return_val=False)

    batch_spinner.succeed("Finished loading from PubChem")
    batch_spinner.stop()


def shred_merge_add_df_mols(dataframe, cmd_pointer):
    """shreds the molecule relevent properties from dataframe and loads into molecules"""

    dict_list = dataframe.to_dict("records")
    merge_list = []
    for i, a_mol in enumerate(dict_list):
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
            merge_mol = get_smol_from_mws(cmd_pointer, name)
            if merge_mol is not None:
                Name_Flag = True
        if not valid_smiles(a_mol["SMILES"]):
            err_msg = f"#{i} - <error>Invalid SMILES, molecule discarded:</error> <yellow>{a_mol['SMILES']}</yellow>"
            output_text(err_msg, return_val=False)
            continue
        merge_mol = get_smol_from_mws(cmd_pointer, a_mol["SMILES"])
        if Name_Flag is True and merge_mol is None:
            # output_error("There is already a molecule by the name, adding  increment to the name  " + name, return_val=False)
            i = 1
            while get_smol_from_mws(cmd_pointer, name + "-" + str(i)) is not None:
                i = i + 1
            name = name + "-" + str(i)

        # if Name_Flag is True and merge_mol["properties"]["canonical_smiles"] != canonical_smiles(a_mol["SMILES"]):
        #    output_error("There is already a molecule by the name, adding  increment to the name " + name, return_val=False)
        #    continue
        if merge_mol is None:
            if name is None:
                name = a_mol["SMILES"]
            merge_mol = new_smol(a_mol["SMILES"], name=name)
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
