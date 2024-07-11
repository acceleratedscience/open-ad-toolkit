"""Store and retrieve results from analysis"""

import os
import glob
import time
import pickle

from openad.molecules.mol_functions import canonicalize
from rdkit import Chem

analysis_record = {"smiles": None, "toolkit": None, "function": None, "parameters": {}, "results": {}}
CACHE_DIR = "/_result_cache/"


def create_analysis_record(smiles, toolkit, function, parameters, results):
    """creates an analysis record for a molecule"""
    record = analysis_record.copy()
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
    except:
        inchikey = result["smiles"]
        return False

    filename = f'{inchikey}-{str(result["toolkit"]).upper()}-{str(result["function"]).upper()}-{timestr}.res'
    _write_analysis(result, _create_workspace_dir_if_nonexistent(cmd_pointer, CACHE_DIR) + filename)
    return True


def _retrieve_results(smiles: str, cmd_pointer) -> list | bool:
    """retrieves results from workspace cache"""
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


def clear_results(cmd_pointer, inp) -> list | bool:
    """clears results from workspace cache"""
    for i in glob.glob(
        _create_workspace_dir_if_nonexistent(cmd_pointer, CACHE_DIR) + "*.res",
        recursive=True,
    ):
        os.remove(i)
    return True


def attach_results(smiles, cmd_pointer) -> bool:
    """attaches a result to an existing molecule in the current working set, if it exists"""

    smiles = canonicalize(smiles)
    results = _retrieve_results(smiles, cmd_pointer)
    for mol in cmd_pointer.molecule_list:
        if "analysis" not in mol:
            mol["analysis"] = []
        if mol["properties"]["canonical_smiles"] == smiles:
            for result in results:
                if canonicalize(mol["properties"]["canonical_smiles"]) == canonicalize(result["smiles"].split("~")[0]):
                    if result not in mol["analysis"]:
                        mol["analysis"].append(result.copy())
    return True


def attach_all_results(cmd_pointer, inp) -> bool:
    """attaches a result to an existing molecule in the current working set, if it exists"""
    i = 0

    while i < len(cmd_pointer.molecule_list):
        mol = cmd_pointer.molecule_list[i].copy()
        if "analysis" not in mol:
            mol["analysis"] = []
        results = _retrieve_results(canonicalize(mol["properties"]["canonical_smiles"]), cmd_pointer)
        for result in results:
            if result not in mol["analysis"] and canonicalize(mol["properties"]["canonical_smiles"]) == canonicalize(
                result["smiles"].split("~")[0]
            ):
                mol["analysis"].append(result)
                cmd_pointer.molecule_list[i] = mol.copy()
        i = i + 1
    return True


def _write_analysis(result: dict, location):
    """writes molecules to a given file"""
    with open(os.path.expanduser(location), "wb") as handle:
        pickle.dump(result, handle)
    return True
