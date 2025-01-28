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
    mws_add,
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

    # Clear mws unless append is passed
    if "append" not in inp:
        cmd_pointer.molecule_list = []

    added_count = 0
    failed_count = 0
    for smol in molset:
        success = mws_add(cmd_pointer, smol, force=True, suppress=True)
        if success:
            added_count += 1
        else:
            failed_count += 1

    # Todo: `load mols using file` should add instead of overwrite your current mols,
    # when this is implemented, we'll need to calculate successfully loaded mols differently.
    if added_count > 0:
        output_success(
            f"Successfully loaded <yellow>{added_count}</yellow> molecules into the working set",
            pad=0,
            return_val=False,
        )
        if failed_count > 0:
            output_error(f"Ignored <yellow>{failed_count}</yellow> duplicates", pad=0, return_val=False)
    else:
        output_error(
            f"No new molecules were added, all {failed_count} provided molecules were are already present in the working set",
            pad=0,
            return_val=False,
        )
    return


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


# # This used to be a part of `load_mols_to_mws` but was replaced by looping mws_add().
# # Ads a result, newly added molecules are not merged, but ignored if they already exist.
# # This merge logic should be added to mws_add() instead. Kept here as a reminder.
# def shred_merge_add_df_mols(molset, cmd_pointer):
#     """Merge properties from a molset into the molecule working set"""
#     # print('shred_merge_add_df_mols')

#     merge_list = []
#     for i, smol in enumerate(molset):

#         identifiers = smol.get("identifiers", {})
#         name = identifiers.get("name")
#         canonical_smiles = identifiers.get("canonical_smiles")
#         # print('\n\n• ', i, canonical_smiles)

#         # Ignore invalid molecules
#         if not valid_smiles(canonical_smiles):
#             output_error(f"#{i} - Invalid SMILES, molecule discarded: <yellow>{smol['SMILES']}</yellow>", return_val=False)
#             continue

#         # Check if molecule with this name already exists in mws
#         name_match = False
#         if name:
#             name_match = bool(get_smol_from_mws(cmd_pointer, name))
#         # print('name_match', name_match)

#         # Check if molecule with this SMILES already exists in mws
#         smiles_match = get_smol_from_mws(cmd_pointer, canonical_smiles)
#         # print('smiles_match', bool(smiles_match))
#         merged_smol = None
#         if smiles_match:
#             merged_smol = merge_smols(smol, smiles_match)
#             print('merged_smol', bool(merged_smol))
#             GLOBAL_SETTINGS["grammar_refresh"] = True

#         # If a molecule with this name already exists in the mws,
#         # but the SMILES don't match, we add an increment to the name.
#         elif name_match:
#             i = 1
#             while get_smol_from_mws(cmd_pointer, f"name-{i}"):
#                 i += 1
#             prev_name = name
#             name = f"name-{i}"
#             output_warning(f"#{i} - Different molecule already named <yellow>{prev_name}</yellow>, renaming to ${name}", return_val=False)
#             continue

#         # print('merged_smol', bool(merged_smol))
#         if merged_smol:
#             j = 0
#             updated_flag = False

#             # Replace the existing molecule with the merged molecule
#             while j < len(cmd_pointer.molecule_list):
#                 mws_smiles = cmd_pointer.molecule_list[j].get("identifiers").get("canonical_smiles")
#                 merged_smol_smiles = smiles_match.get("identifiers").get("canonical_smiles")
#                 if (merged_smol_smiles == mws_smiles):
#                     cmd_pointer.molecule_list[j] = merged_smol
#                     merged_smol_name = smiles_match.get("identifiers").get("name")
#                     if not merged_smol_name:
#                         merged_smol_name = merged_smol_smiles
#                     merge_list.append(merged_smol_name)
#                     j = len(cmd_pointer.molecule_list)
#                     updated_flag = True
#                 j += 1

#             # If the molecule was not found in the mws, add it
#             if updated_flag is False:
#                 cmd_pointer.molecule_list.append(smiles_match)

#     if len(merge_list) > 0:
#         merge_list = list(set(merge_list))
#         merge_list.sort()
#         output_warning("The following molecules had one or more duplicate entries and were merged:", return_val=False)
#         print(merge_list)
#         # output_text("\n- " + "\n- ".join(merge_list), return_val=False)
