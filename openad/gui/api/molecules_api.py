import os
import json
import time
import shutil
import asyncio
import aiofiles
from rdkit import Chem
from flask import Response, request
from urllib.parse import unquote

# from openad.molecules.mol_api import get_molecule_data
from openad.molecules.mol_commands import retrieve_mol
from openad.molecules.mol_functions import organize_properties
from openad.molecules.mol_functions import mol2svg, mol2sdf

from openad.helpers.json_decimal_encoder import DecimalEncoder
from openad.helpers.files import open_file
from openad.helpers.output import output_error


class MoleculesApi:
    def __init__(self, cmd_pointer):
        self.cmd_pointer = cmd_pointer

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

        # Success
        if mol:
            mol = organize_properties(mol)
            return mol

        # Fail
        else:
            response = Response(None, status=500)
            response.status = f"No molecule found with provided identifier '{identifier}'"
            return response

    # Get a molecule's SDF and SVG.
    # Used when opening a .mol.json file.
    def get_mol_viz_data(self):
        data = json.loads(request.data) if request.data else {}
        inchi_or_smiles = data["inchi_or_smiles"]

        mol_rdkit = Chem.MolFromInchi(inchi_or_smiles)
        if not mol_rdkit:
            mol_rdkit = Chem.MolFromSmiles(inchi_or_smiles)  # pylint: disable=no-member (false positive)
        if mol_rdkit:
            svg = mol2svg(mol_rdkit)
            sdf = mol2sdf(mol_rdkit)
        else:
            svg, sdf = None, None

        return {"sdf": sdf, "svg": svg}

    # Filter a molecule set by a string.
    def get_molset(self):
        print("get_molset")
        data = json.loads(request.data) if request.data else {}
        path = unquote(data["path"]) if "path" in data else ""
        cache_id = data["cacheId"] if "cacheId" in data else ""
        query = data["query"] if "query" in data else {}

        # Parse the query
        search_str = query["search"] if "search" in query else ""
        page = int(query["page"]) if "page" in query else 1
        sort = query["sort"] if "sort" in query else None
        page_size = 48  # Hardcoded for now

        # When we open a molset, it is "fresh".
        # We have to assign it a cash_id and store a working copy.
        fresh = not cache_id

        # Workspace path
        workspace_path = self.cmd_pointer.workspace_path(self.cmd_pointer.settings["workspace"])

        # Opening a file: assign a cache_id and store a working copy.
        if fresh:
            cache_id = int(time.time() * 1000)
            file_path = workspace_path + "/" + path

            # Make a working copy of the molset.
            if os.path.exists(file_path):
                cache_path = self._cache_path(cache_id)
                shutil.copy(file_path, cache_path)

                # Add indices to molecules in our working copy,
                # without blocking the thread.
                asyncio.run(self._add_index(cache_path))

        else:
            cache_path = self._cache_path(cache_id)

        # Read file from cache.
        molset, err_code = open_file(cache_path, return_err="code", dumb=False)  # %%

        # Return error
        if err_code:
            return err_code, 500

        # Only add index when the working copy hasn't been indexed yet.
        if fresh:
            print("@@ ADD")
            for i, mol in enumerate(molset):
                mol["index"] = i + 1

        # Search all identifiers and properties for the query.
        if search_str:
            results = []
            for mol in molset:
                found = False
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
            if (sort) == "name":
                results = sorted(results, key=lambda mol: self._sort_mol(mol, sort, "identifiers"))
            elif sort:
                results = sorted(results, key=lambda mol: self._sort_mol(mol, sort, "properties"))
        # In the edge case where our dataset mixes string and
        # number values, we want to avoid crashing the app.
        except TypeError as err:
            output_error(err)

        # Store matching ids before pagination
        matching = [mol["index"] for mol in results]

        # Paginate
        total_pages = len(results) // page_size + 1
        page = min(page, total_pages)  # Make sure that page number is lowered in case of a too high value.
        total = len(results)
        skip = 48 * (page - 1)
        results = results[skip : skip + 48]
        # result_page = json.dumps(results[skip : skip + 48], cls=DecimalEncoder) # Stringified, not needed

        return {
            "cacheId": cache_id,  # Used to identify our working copy in next requests
            "mols": results,  # One page of molecules
            "matching": matching,  # Ids of ALL matching molecules
            # Query parameters:
            "searchStr": search_str,
            "sort": sort,
            "total": total,
            "page": page,
            "pageSize": page_size,
        }

    # Add an index to our cached working copy, without blocking the thread.
    # (This can take up to a second for very large files, say ~100MB)
    # Note: read and write happen separately to make sure the file is overwritten.
    async def _add_index(self, cache_path):
        start_time = time.time()
        # Read
        async with aiofiles.open(cache_path, "r", encoding="utf-8") as f:
            content = await f.read()
        molset = json.loads(content)
        for i, mol in enumerate(molset):
            mol["index"] = i + 1
        # Write
        async with aiofiles.open(cache_path, "w", encoding="utf-8") as f:
            await f.write(json.dumps(molset, ensure_ascii=False, indent=4, cls=DecimalEncoder))
        print("_add_index took %s seconds" % (time.time() - start_time))

    # Sorter function.
    def _sort_mol(self, mol, sort, category):
        value = (mol.get(category) or {}).get(sort)
        # Returning a tuple will sort by the first value, then the second, etc.
        # This lets us group all none values on top.
        return (value is None, value.lower() if isinstance(value, str) else value)

    # Remove molecules from a molset cache
    def remove_from_molset(self):
        print("remove_from_molset")
        data = json.loads(request.data) if request.data else {}
        cache_id = data["cacheId"] if "cacheId" in data else ""
        indices = data["indices"] if "indices" in data else []

        if not cache_id or not indices:
            return f"Unrecognized cache_id: {cache_id}", 500

        # Compile path
        cache_path = self._cache_path(cache_id)

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

        print("@@@")
        return self.get_molset()

    # Clear a molset's cached working copy.
    def clear_from_cache(self):
        data = json.loads(request.data) if request.data else {}
        cache_id = data["cacheId"] if "cacheId" in data else ""
        print("clear_from_cache", data, cache_id)

        if not cache_id:
            return f"Unrecognized cache_id: {cache_id}", 500

        cache_path = self._cache_path(cache_id)
        if os.path.exists(cache_path):
            os.remove(cache_path)

        return "ok"

    # Compile the file path to the cached working copy of a molset.
    def _cache_path(self, cache_id):
        return (
            self.cmd_pointer.workspace_path(self.cmd_pointer.settings["workspace"])
            + "/.cache/molset-"
            + str(cache_id)
            + ".json"
        )
