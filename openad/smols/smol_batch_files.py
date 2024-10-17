""" handles back load and unloading for molecule operations"""

import pandas
from rdkit import RDLogger
from rdkit.Chem import PandasTools
from openad.smols.smol_functions import (
    get_smol_from_mws,
    merge_molecule_properties,
    valid_smiles,
    get_smol_from_pubchem,
    new_smol,
    canonicalize,
    load_mols_from_file,
    merge_smols,
    mws_remove,
)
from openad.smols.smol_transformers import dataframe2molset
from openad.app.global_var_lib import GLOBAL_SETTINGS
from openad.helpers.output import output_error, output_warning, output_success, output_text
from openad.helpers.output_msgs import msg
from openad.helpers.spinner import Spinner
from openad.plugins.style_parser import style

# Suppress RDKit errors
RDLogger.DisableLog("rdApp.error")
RDLogger.DisableLog("rdApp.warning")

mol_name_cache = {}


def merge_molecule_property_data(cmd_pointer, inp=None, dataframe=None):
    """
    Merge data from a dataframe into your molecule working set molecule properties.

    The dataframe should contain the following columns:
    - SMILES/subject
    - property
    - value
    """

    if dataframe is None and inp is None:
        return False

    if dataframe is None:
        # Load from dataframe
        if (
            "merge_molecules_data_dataframe" in inp.as_dict()
            or "merge_molecules_data_dataframe-DEPRECATED"  # Can be removed once the deprecated syntax has been removed
            in inp.as_dict()
        ):
            dataframe = cmd_pointer.api_variables[inp.as_dict()["in_dataframe"]]

        # Load from file (not yet implemented)
        else:
            dataframe = load_mol_data(inp.as_dict()["moles_file"], cmd_pointer)

        if dataframe is None:
            output_error("Source not found ", return_val=False)
            return True

    # Detect the SMILES/subject column
    if "subject" in dataframe.columns:
        smiles_key = "subject"
    elif "smiles" in dataframe.columns:
        smiles_key = "smiles"
    elif "SMILES" in dataframe.columns:
        smiles_key = "SMILES"
    else:
        output_error(
            "No <yellow>subject</yellow> or <yellow>SMILES</yellow> column found in merge data", return_val=False
        )
        return True

    # Detect the property column
    if "property" in dataframe.columns:
        prop_key = "property"
    elif "PROPERTY" in dataframe.columns:
        prop_key = "PROPERTY"
    else:
        output_error("No <yellow>property</yellow> column found in merge data", return_val=False)
        return True

    # Detect the value column
    if "value" in dataframe.columns:
        val_key = "value"
    elif "VALUE" in dataframe.columns:
        val_key = "VALUE"
    elif "result" in dataframe.columns:
        val_key = "result"
    elif "RESULT" in dataframe.columns:
        val_key = "RESULT"
    else:
        output_error("No <yellow>result</yellow> or <yellow>value</yellow> column found", return_val=False)
        return True

    # Pivot the dataframe
    dataframe = dataframe.pivot_table(index=smiles_key, columns=[prop_key], values=val_key, aggfunc="first")
    dataframe = dataframe.reset_index()

    for row in dataframe.to_dict("records"):
        update_flag = True
        merge_smol = None

        try:
            smiles = canonicalize(row[smiles_key])
            merge_smol = get_smol_from_mws(cmd_pointer, smiles)
            GLOBAL_SETTINGS["grammar_refresh"] = True
        except Exception:  # pylint: disable=broad-except
            output_warning("unable to canonicalise:" + row[smiles_key])
            continue

        if merge_smol is None:
            merge_smol = new_smol(smiles, name=row[smiles_key])
            update_flag = False
        else:
            update_flag = True

        if merge_smol is not None:
            smol = merge_molecule_properties(row, merge_smol)
            GLOBAL_SETTINGS["grammar_refresh"] = True
            if update_flag is True:
                mws_remove(cmd_pointer, merge_smol, force=True, suppress=True)
            cmd_pointer.molecule_list.append(smol)

    output_success("Data merged into your working set", return_val=False)
    GLOBAL_SETTINGS["grammar_refresh"] = True
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


# Unused
# def load_molset_to_mws(cmd_pointer: object, molset: list, append=False):
#     """
#     Load a molset into the molecule working set.
#     """
#     if append:
#         cmd_pointer.molecule_list.extend(molset)
#     else:
#         cmd_pointer.molecule_list = molset


def load_mols_to_mws(cmd_pointer, inp):
    """
    Load a batch of molecules into the molecule working set.
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

    # Add PubChem data
    if "enrich_pubchem" in inp.as_dict():
        _enrich_with_pubchem_data(cmd_pointer, molset)

    # Append or overwrite molecules
    if "append" in inp:
        cmd_pointer.molecule_list.extend(molset)
    else:
        cmd_pointer.molecule_list = molset

    if molset:
        # shred_merge_add_df_mols(molset, cmd_pointer)
        # Todo: `load mols using file` should add instead of overwrite your current mols,
        # when this is implemented, we'll need to calculate successfully loaded mols differently.
        return output_success(
            f"Successfully loaded <yellow>{len(molset)}</yellow> molecules into the working set", pad=0
        )


def _enrich_with_pubchem_data(cmd_pointer, molset):
    """
    Pull data from PubChem to merge in into a molset.
    """

    output_molset = []

    spinner = Spinner(GLOBAL_SETTINGS["VERBOSE"])
    spinner.start("Fetching from PubChem")

    for i, smol in enumerate(molset):
        try:
            identifiers = smol["identifiers"]

            # Get name field regardless of case
            name = next((value for key, value in identifiers.items() if key.lower() == "name"), None)
            spinner.text = style(f"<soft>Fetching from PubChem: #{i} - {name}</soft>")

            # Use fallback name is missing
            if not name:
                name = smol.get("chemical_name", None)

            # Select the identifier keys we'll look for in order of preference
            keys = ["inchi", "canonical_smiles", "isomeric_smiles", "smiles", "inchikey", "name", "cid"]
            identifier = next((identifiers.get(key) for key in keys if identifiers.get(key) is not None), None)
            name = name or identifier or "Unknown molecule"
            if not identifier:
                output_warning(f"#{i} - No valid identifier found for {name}", return_val=False)
                continue

            # Fetch enriched molecule
            smol_enriched = get_smol_from_pubchem(identifier)
            if not smol_enriched:
                output_warning(f"#{i} - Failed to enrich {name}", return_val=False)

            # Merge enriched data
            smol = merge_smols(smol, smol_enriched)
            output_molset.append(smol)

        except Exception as err:  # pylint: disable=broad-except
            spinner.stop()
            output_error(["Something went wrong enriching molecules with data from PubChem", err], return_val=False)

    spinner.succeed("Done")
    spinner.stop()
    return output_molset


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

        # if Name_Flag is True and merge_mol["identifiers"]["canonical_smiles"] != canonical_smiles(a_mol["SMILES"]):
        #    output_error("There is already a molecule by the name, adding  increment to the name " + name, return_val=False)
        #    continue
        if merge_mol is None:
            if name is None:
                name = a_mol["SMILES"]
            merge_mol = new_smol(a_mol["SMILES"], name=name)
        if merge_mol is not None:
            merge_mol = merge_molecule_properties(a_mol, merge_mol)
            GLOBAL_SETTINGS["grammar_refresh"] = True
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
                merge_mol["identifiers"]["canonical_smiles"]
                == cmd_pointer.molecule_list[i]["identifiers"]["canonical_smiles"]
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
