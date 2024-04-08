import os
import json
import time
import copy
import shutil
import asyncio
import aiofiles
from rdkit import Chem
from rdkit.Chem import PandasTools
from flask import Response, request
from urllib.parse import unquote

# from openad.molecules.mol_api import get_molecule_data
from openad.molecules.mol_commands import retrieve_mol
from openad.molecules.mol_functions import new_molecule, molformat_v2, mol2svg, mol2sdf

# from openad.workers.file_system import fs_assemble_cache_path
import openad.workers.file_system as fs

from openad.helpers.json_decimal_encoder import DecimalEncoder
from openad.helpers.files import open_file
from openad.helpers.output import output_error


class MoleculesApi:
    def __init__(self, cmd_pointer):
        self.cmd_pointer = cmd_pointer

    # -----------------------------
    # Molecules
    # -----------------------------

    # Get molecule data, plus SDF and SVG.
    # Used when requesting a molecule by its identifier.
    def get_mol(self):
        data = json.loads(request.data) if request.data else {}
        identifier = data["identifier"]

        if not identifier:
            response = Response(None, status=500)
            response.status = "No identifier provided."
            return response

        mol = retrieve_mol(identifier)

        # Fail
        if not mol:
            response = Response(None, status=500)
            response.status = f"No molecule found with provided identifier '{identifier}'"
            return response

        # Success
        else:
            mol = molformat_v2(mol)
            return mol, 200

    # Get a molecule's SDF and SVG.
    # Used when opening a .mol.json file.
    def get_mol_viz_data(self):
        data = json.loads(request.data) if request.data else {}
        inchi_or_smiles = data["inchi_or_smiles"] if "inchi_or_smiles" in data else None

        if not inchi_or_smiles:
            return "get_mol_viz_data() -> Invalid inchi_or_smiles", 500

        mol_rdkit = Chem.MolFromInchi(inchi_or_smiles)
        if not mol_rdkit:
            mol_rdkit = Chem.MolFromSmiles(inchi_or_smiles)  # pylint: disable=no-member (false positive)
        if mol_rdkit:
            svg = mol2svg(mol_rdkit)
            sdf = mol2sdf(mol_rdkit)
        else:
            svg, sdf = None, None

        return {"sdf": sdf, "svg": svg}, 200

    def get_mol_data_from_molset(self):
        data = json.loads(request.data) if request.data else {}
        cache_id = data["cacheId"] if "cacheId" in data else ""
        index = data["index"] if "index" in data else 1

        # Workspace path
        cache_path = fs.fs_assemble_cache_path(self.cmd_pointer, "molset", cache_id)

        # Read file from cache.
        molset, err_code = open_file(cache_path, return_err="code")

        # Return error
        if err_code:
            return err_code, 500

        return molset[index - 1], 200

    # -----------------------------
    # Molsets
    # -----------------------------

    # Filter a molecule set by a string.
    def query_molset(self):
        # print("get_molset")
        data = json.loads(request.data) if request.data else {}
        cache_id = data["cacheId"] if "cacheId" in data else ""
        query = data["query"] if "query" in data else {}

        # Read file from cache.
        cache_path = fs.fs_assemble_cache_path(self.cmd_pointer, "molset", cache_id)
        molset, err_code = get_molset_mols(cache_path)

        # Return error
        if err_code:
            return err_code, 500

        return create_molset_response(molset, query, cache_id), 200

    # Get your working list of molecules.
    def get_my_mols(self):
        data = json.loads(request.data) if request.data else {}
        query = data["query"] if "query" in data else {}

        if len(self.cmd_pointer.molecule_list) > 0:
            molset = []
            for i, mol in enumerate(self.cmd_pointer.molecule_list):
                mol_v2 = molformat_v2(mol)
                mol_v2["index"] = i + 1
                molset.append(mol_v2)
            return create_molset_response(molset, query), 200
        else:
            return "No molecules in list", 204

    # Remove molecules from a molset cache
    def remove_from_molset(self):
        # print("remove_from_molset")
        data = json.loads(request.data) if request.data else {}
        cache_id = data["cacheId"] if "cacheId" in data else ""
        indices = data["indices"] if "indices" in data else []
        query = data["query"] if "query" in data else {}

        if not cache_id:
            return f"remove_from_molset() -> Unrecognized cache_id: {cache_id}", 500

        if len(indices) == 0:
            return f"remove_from_molset() -> No indices provided", 500

        # Compile path
        cache_path = fs.fs_assemble_cache_path(self.cmd_pointer, "molset", cache_id)

        # Read file from cache
        molset, err_code = open_file(cache_path, return_err="code", dumb=False)  # %%

        # Return error
        if err_code:
            return err_code, 500

        # Remove molecules
        molset = [mol for mol in molset if mol["index"] not in indices]

        # Write to cache
        with open(cache_path, "w", encoding="utf-8") as f:
            json.dump(molset, f, ensure_ascii=False, indent=4, cls=DecimalEncoder)

        # Create response object.
        return create_molset_response(molset, query, cache_id), 200

    def save_molset_changes(self):
        # print("save_molset_changes")
        data = json.loads(request.data) if request.data else {}
        cache_id = data["cacheId"] if "cacheId" in data else ""
        path = unquote(data["path"]) if "path" in data else ""

        # Compile path
        workspace_path = self.cmd_pointer.workspace_path(self.cmd_pointer.settings["workspace"])
        file_path = workspace_path + "/" + path

        if not cache_id:
            return f"save_molset_changes() -> Unrecognized cache_id: {cache_id}", 500

        cache_path = fs.fs_assemble_cache_path(self.cmd_pointer, "molset", cache_id)

        if not os.path.exists(file_path):
            return f"save_molset_changes() -> File not found: {file_path}", 500

        if not os.path.exists(cache_path):
            return f"save_molset_changes() -> Cached working copy not found: {cache_path}", 500

        shutil.copy(cache_path, file_path)
        # os.remove(cache_path)

        return "ok", 200

    # Clear a molset's cached working copy.
    def clear_molset_working_copy(self):
        # print("clear_molset_working_copy", data, cache_id)
        data = json.loads(request.data) if request.data else {}
        cache_id = data["cacheId"] if "cacheId" in data else ""

        if not cache_id:
            return f"clear_molset_working_copy() -> Unrecognized cache_id: {cache_id}", 500

        cache_path = fs.fs_assemble_cache_path(self.cmd_pointer, "molset", cache_id)
        if os.path.exists(cache_path):
            os.remove(cache_path)

        return "ok", 200


# Return molset from JSON file
def get_molset_mols(path_absolute):
    """
    Return the list of molecules from a molset file,
    with an index added to each molecule.

    Parameters
    ----------
    path_absolute: str
        The absolute path to the molset file.

    Returns: data, err_code
    """
    # Read file contents
    molset, err_code = open_file(path_absolute, return_err="code")

    # Add index.
    if molset:
        for i, mol in enumerate(molset):
            mol["index"] = i + 1

    return molset, err_code


# Return molset from SMILES file
def smiles_file2molset(path_absolute):
    """
    Read the content of a .smi file and return a molset.
    Specs for .smi files: http://opensmiles.org/opensmiles.html - 4.5

    This takes about 3 seconds per 10,000 molecules on an Apple M2 with 16GB or memory.
    """
    # Read file's content
    data, err_code = open_file(path_absolute, return_err="code")
    if err_code:
        return None, err_code

    # Parse SMILES
    smiles_list = data.splitlines()
    # Ignore any properties that may be listed after the SMILES string.
    smiles_list = [smiles.split(" ")[0] for smiles in smiles_list if smiles]
    molset = []
    for i, smiles in enumerate(smiles_list):
        mol = new_molecule(smiles)
        if mol:
            mol = molformat_v2(mol)
        else:
            mol = {
                "identifiers": {"canonical_smiles": smiles},
                "properties": {},
            }

        mol["index"] = i + 1
        molset.append(mol)

    return molset, None


def sdf_path2molset(sdf):
    from openad.molecules.mol_functions import OPENAD_MOL_DICT

    try:
        mols = Chem.SDMolSupplier(sdf)  # pylint: disable=no-member
        molset = []
        for i, mol in enumerate(mols):
            mol_dict = copy.deepcopy(OPENAD_MOL_DICT)
            mol_dict["properties"] = {prop: mol.GetProp(prop) for prop in mol.GetPropNames()}
            mol_dict = molformat_v2(mol_dict)
            mol_dict["index"] = i + 1
            molset.append(mol_dict)
        return molset, None
    except Exception as err:
        return None, err


def sdf_data2molset(sdf_data):
    from openad.molecules.mol_functions import OPENAD_MOL_DICT

    try:
        # Split the SDF data into blocks and convert each block to a Mol object
        sdf_blocks = sdf_data.split("\n$$$$\n")
        mols = [Chem.MolFromMolBlock(block) for block in sdf_blocks if block]  # pylint: disable=no-member

        molset = []
        for i, mol in enumerate(mols):
            mol_dict = copy.deepcopy(OPENAD_MOL_DICT)
            mol_dict["properties"] = {prop: mol.GetProp(prop) for prop in mol.GetPropNames()}
            mol_dict = molformat_v2(mol_dict)
            mol_dict["index"] = i + 1
            molset.append(mol_dict)
        return molset, None
    except Exception as err:
        return None, err


def df2sdf(df):
    """
    Reads a dataframe, looks for an InChI or
    SMILES column and returns a molset.
    """
    # Create an empty DataFrame with an 'ROMol' column
    # PandasTools.AddMoleculeColumnToFrame(df, smilesCol="SMILES", molCol="ROMol")

    colsLowercase = [col.lower() for col in df.columns]

    if "inchi" in colsLowercase:
        index = colsLowercase.index("inchi")
        key = df.columns[index]
        key_type = "inchi"
    elif "smiles" in colsLowercase:
        index = colsLowercase.index("smiles")
        key = df.columns[index]
        key_type = "smiles"
    else:
        return None

    # Convert the molecules to SDF format
    sdf_data = ""
    for i, row in df.iterrows():
        if key_type == "inchi":
            mol_rdkit = Chem.MolFromInchi(row[key])  # pylint: disable=no-member
        elif key_type == "smiles":
            mol_rdkit = Chem.MolFromSmiles(row[key])  # pylint: disable=no-member

        if mol_rdkit is not None:
            sdf_data += Chem.MolToMolBlock(mol_rdkit) + "\n$$$$\n"  # pylint: disable=no-member

    return sdf_data


def df2molset(df):
    sdf_data = df2sdf(df)
    if sdf_data is None:
        return None
    return sdf_data2molset(sdf_data)


def create_molset_response(molset, query={}, cache_id=None):
    """
    Return a filtered and paginated subset of a molset, wrapped into
    a response object that is ready to be consumed by the frontend.

    Parameters
    ----------
    molset: list
        The entire list of molecules from a file.
    cache_id: str
        The unique identifier of the working copy of the molset.
    query: dict
        The query parameters from the frontend.
        Eg. { queryStr: "C1=CC=CC=C1", smarts: 1, page: 1, sort: "name" }
    """

    # Parse the query
    search_str = unquote(query["search"]) if "search" in query else ""
    smarts_mode = query["smarts"] if "smarts" in query else False
    page = int(query["page"]) if "page" in query else 1
    sort = query["sort"] if "sort" in query else None
    page_size = 48  # Hardcoded for now

    # Filter by query
    if search_str:
        results = []
        for mol in molset:
            found = False

            # Substructure search - match against smiles only.
            if smarts_mode:
                if search_str.lower() in mol["identifiers"]["canonical_smiles"].lower():
                    results.append(mol)

            # Regular search - match against all identifiers and properties.
            else:
                for key in mol["identifiers"]:
                    if search_str.lower() in str(mol["identifiers"][key]).lower():
                        results.append(mol)
                        found = True
                        break
                if not found:
                    for key in mol["properties"]:
                        if search_str.lower() in str(mol["properties"][key]).lower():
                            results.append(mol)
                            break

    # No query
    else:
        results = molset

    # Sort
    try:
        # print("INDEX:", molset[0].get("index"))
        # print("SORT:", sort)
        if sort:
            reverse = sort.startswith("-")
            sort_key = sort.lstrip("-")
            results = sorted(results, key=lambda mol: _sort_mol(mol, sort_key), reverse=reverse)
    except TypeError as err:
        # In the edge case where our dataset mixes string and
        # number values, we want to avoid crashing the app.
        output_error(err)

    # Store matching ids before pagination.
    matching = [mol.get("index") for mol in results]

    # Paginate
    total_pages = len(results) // page_size + 1
    page = min(page, total_pages)  # Make sure that page number is lowered in case of a too high value.
    total = len(molset)
    result_count = len(results)
    skip = 48 * (page - 1)
    results = results[skip : skip + 48]

    return {
        "cacheId": cache_id,  # Used to identify our working copy in next requests
        "mols": results,  # One page of molecules
        "matching": matching,  # Ids of ALL matching molecules
        "total": total,
        "resultCount": result_count,
        # Query parameters:
        "searchStr": search_str,
        "searchMode": "smarts" if smarts_mode else "text",
        "sort": sort,
        "page": page,
        "pageSize": page_size,
    }


def _sort_mol(mol, sort_key):
    """
    Sorter function for a molset.

    Parameters
    ----------
    mol: dict
        A molecule object.
    sort_key: str
        The key of the category whose value we'll sort by.
        Eg. 'name' (identifier) or 'molecular_weight' (property).
    """
    if sort_key == "index":
        value = mol.get(sort_key)
    elif sort_key == "name":
        value = (mol.get("identifiers") or {}).get(sort_key)
    else:
        value = (mol.get("properties") or {}).get(sort_key)

    value = __prep_sort_value(value)

    # Returning a tuple will sort by the first value, then the
    # second, etc. This lets us group all none values on top.
    return (value is None, value)


def __prep_sort_value(value):
    """
    Prepare a value for sorting by:
    - Converting strings to lowercase
    - Converting number strings to floats

    Parameters
    ----------
    value: str, int, float
        The value to prepare.
    """

    # Convert number strings to floats
    try:
        return float(value)
    except Exception as err:
        pass

    # Convert strings to lowercase
    if isinstance(value, str):
        return value.lower()

    return value


def index_molset_file_async(path_absolute):
    """
    Add an index to each molecule of a molset file,
    without blocking the main thread.

    This is used to index a cached working copy of a molset
    right after it's created.

    Parameters
    ----------
    cache_path: str
        The path to the cached working copy of a molset.
    """

    async def _index_molset_file(cache_path):
        # Read
        async with aiofiles.open(cache_path, "r", encoding="utf-8") as f:
            content = await f.read()
        molset = json.loads(content)
        for i, mol in enumerate(molset):
            mol["index"] = i + 1
        # Write
        async with aiofiles.open(cache_path, "w", encoding="utf-8") as f:
            await f.write(json.dumps(molset, ensure_ascii=False, indent=4, cls=DecimalEncoder))

    asyncio.run(_index_molset_file(path_absolute))
