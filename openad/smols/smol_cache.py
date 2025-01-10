"""Store and retrieve results from analysis"""

import os
import glob
import copy
import time
import pickle
from rdkit import Chem
from openad.smols.smol_functions import canonicalize
from openad.helpers.data_formats import ANALYSIS_RECORD
from openad.helpers.output import output_success, output_error, output_warning


CACHE_DIR = "/_result_cache/"


# To do: rename toolkit param to plugin (currently consumed by both)
def create_analysis_record(smiles, toolkit, function, parameters, results):
    """creates an analysis record for a molecule"""
    record = copy.deepcopy(ANALYSIS_RECORD)
    record["smiles"] = smiles
    record["toolkit"] = toolkit
    record["function"] = function
    record["parameters"] = parameters
    record["results"] = results
    return record


def _create_workspace_dir_if_nonexistent(cmd_pointer, dir_name):
    if not os.path.isdir(cmd_pointer.workspace_path(cmd_pointer.settings["workspace"].upper()) + "/" + dir_name):
        os.mkdir(cmd_pointer.workspace_path(cmd_pointer.settings["workspace"].upper()) + "/" + dir_name)
    return cmd_pointer.workspace_path(cmd_pointer.settings["workspace"].upper()) + "/" + dir_name


def save_result(result: dict, cmd_pointer) -> bool:
    """saves a result of an analysis tool"""

    try:
        timestr = time.strftime("%Y%m%d-%H%M%S")

        rdkit_mol = Chem.MolFromSmiles(result["smiles"])  # pylint: disable=no-member

        inchi = Chem.rdinchi.MolToInchi(rdkit_mol)[0]

        inchikey = Chem.inchi.InchiToInchiKey(inchi)
    except Exception:  # pylint: disable=broad-except
        inchikey = result["smiles"]
        return False

    filename = f'{inchikey}-{str(result["toolkit"]).upper()}-{str(result["function"]).upper()}-{timestr}.res'
    _write_analysis(result, _create_workspace_dir_if_nonexistent(cmd_pointer, CACHE_DIR) + filename)
    return True


def _retrieve_results(smiles: str, cmd_pointer) -> list | bool:
    """
    Retrieve results from workspace cache
    """

    rdkit_mol = Chem.MolFromSmiles(smiles)  # pylint: disable=no-member
    inchi = Chem.rdinchi.MolToInchi(rdkit_mol)[0]
    inchikey = Chem.inchi.InchiToInchiKey(inchi)
    results = []
    for i in glob.glob(
        _create_workspace_dir_if_nonexistent(cmd_pointer, CACHE_DIR) + inchikey + "-*.res",
        recursive=True,
    ):
        func_file = open(i, "rb")
        results.append(pickle.load(func_file))

    return results.copy()


def results_exist(cmd_pointer) -> bool:
    """
    Check if any results exist in the workspace cache.
    """

    if len(glob.glob(_create_workspace_dir_if_nonexistent(cmd_pointer, CACHE_DIR) + "*.res", recursive=True)) > 0:
        return True
    return False


def clear_analysis(cmd_pointer, inp) -> list | bool:
    """clears results from workspace cache"""

    # No alalysis results
    if results_exist(cmd_pointer) is False:
        return output_warning("Analysis results were already empty")

    # Clear
    for i in glob.glob(
        _create_workspace_dir_if_nonexistent(cmd_pointer, CACHE_DIR) + "*.res",
        recursive=True,
    ):
        os.remove(i)

    # Success
    return output_success("Analysis results have been cleared", return_val=True)


def attach_results(smiles, cmd_pointer) -> bool:
    """attaches a result to an existing molecule in the current working set, if it exists"""

    smiles = canonicalize(smiles)
    results = _retrieve_results(smiles, cmd_pointer)
    for mol in cmd_pointer.molecule_list:
        if "analysis" not in mol:
            mol["analysis"] = []
        if mol["identifiers"]["canonical_smiles"] == smiles:
            for result in results:
                if canonicalize(mol["identifiers"]["canonical_smiles"]) == canonicalize(result["smiles"].split("~")[0]):
                    if result not in mol["analysis"]:
                        mol["analysis"].append(result.copy())
    return True


def enrich_mws_with_analysis(cmd_pointer, inp) -> bool:
    """
    Attach the latest analysis results to a molecule in the current working set, if it present.
    """

    mws = cmd_pointer.molecule_list
    enriched_smols_count = 0

    # No molecules
    if len(mws) == 0:
        return output_error("No molecules in the working set")

    # No alalysis results
    if results_exist(cmd_pointer) is False:
        return output_error("No analysis results found in the workspace cache")

    i = 0
    while i < len(mws):
        smol = mws[i].copy()
        if "analysis" not in smol:
            smol["analysis"] = []
        results = _retrieve_results(canonicalize(smol["identifiers"]["canonical_smiles"]), cmd_pointer)
        for result in results:
            if result not in smol["analysis"] and canonicalize(smol["identifiers"]["canonical_smiles"]) == canonicalize(
                result["smiles"].split("~")[0]
            ):
                smol["analysis"].append(result)
                mws[i] = smol.copy()
                enriched_smols_count = enriched_smols_count + 1
        i = i + 1

    # Success
    if enriched_smols_count > 0:
        return output_success(
            [
                f"{enriched_smols_count}/{len(mws)} molecules in your working set have been enriched with the latest analysis results",
                "Run <cmd>show mols</cmd> to view the updated working set",
            ]
        )

    # Error
    else:
        return output_error("No matching analysis results found for any of the molecules in your working set")


def _write_analysis(result: dict, location):
    """writes molecules to a given file"""
    with open(os.path.expanduser(location), "wb") as handle:
        pickle.dump(result, handle)
    return True
