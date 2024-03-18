import os
import json
import time
import shutil
from rdkit import Chem
from flask import Response, request
from urllib.parse import unquote

# from openad.molecules.mol_api import get_molecule_data
from openad.molecules.mol_commands import retrieve_mol
from openad.molecules.mol_functions import organize_properties
from openad.molecules.mol_functions import mol2svg, mol2sdf

# from openad.helpers.json_decimal_encoder import DecimalEncoder
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
        data = json.loads(request.data) if request.data else {}
        path = unquote(data["path"]) if "path" in data else ""
        cache_id = data["cacheId"] if "cacheId" in data else ""
        query = data["query"] if "query" in data else ""
        page = data["page"] if "page" in data else 1
        page_size = data["pageSize"] if "pageSize" in data else 48
        sort = data["sort"] if "sort" in data else None
        index_only = data["indexOnly"] if "indexOnly" in data else False

        # Workspace path
        workspace_path = self.cmd_pointer.workspace_path(self.cmd_pointer.settings["workspace"])
        filename = path.split("/")[-1].replace(".molset.json", "")

        # Opening a file: store cached copy and assign cache_id.
        if not cache_id:
            cache_id = int(time.time() * 1000)
            print(cache_id)
            file_path = workspace_path + "/" + path
            if os.path.exists(file_path):
                cache_path = workspace_path + "/.cache/" + filename + "--" + str(cache_id) + "molset.json"
                shutil.copy(file_path, cache_path)
        else:
            cache_path = workspace_path + "/.cache/" + filename + "--" + str(cache_id) + "molset.json"

        # Read file from cache.
        molset, err_code = open_file(cache_path, return_err="code", dumb=False)  # %%

        # Return error
        if err_code:
            return None

        # Search all identifiers and properties for the query.
        if query:
            results = []
            for i, mol in enumerate(molset):
                mol["index"] = i + 1
                found = False
                for key in mol["identifiers"]:
                    if query.lower() in str(mol["identifiers"][key]).lower():
                        results.append(mol)
                        break
                if found:
                    continue
                for key in mol["properties"]:
                    if query.lower() in str(mol["properties"][key]).lower():
                        results.append(mol)
                        break

        # No query
        else:
            for i, mol in enumerate(molset):
                mol["index"] = i + 1
            results = molset

        # Sort
        try:
            if (sort) == "name":
                results = sorted(results, key=lambda mol: _sort_mol(mol, sort, "identifiers"))
            elif sort:
                results = sorted(results, key=lambda mol: _sort_mol(mol, sort, "properties"))
        except TypeError as err:
            # In the edge case our dataset mixes string and
            # number values, we want to avoid crashing the app.
            output_error(err)

        # Store matching ids before pagination
        matching = [mol["index"] for mol in results]

        # Paginate
        if not index_only:
            total = len(results)
            skip = 48 * (page - 1)
            # result_page = json.dumps(results[skip : skip + 48], cls=DecimalEncoder) # Stringified, not needed
            results = results[skip : skip + 48]

        # Result indeces only.
        # This is used to select matching molecules past current page.
        if index_only:
            data = {
                "indices": [mol["index"] for mol in results],
            }

        # Regular results.
        else:
            data = {
                "cacheId": cache_id,
                "mols": results,
                "matching": matching,
                "total": total,
                "page": page,
                "pageSize": page_size,
            }

        return data

    # Remove molecules from a molset cache.
    def remove_from_molset(self):
        data = json.loads(request.data) if request.data else {}
        cache_id = data["cacheId"] if "cacheId" in data else ""
        indices = data["indices"] if "indices" in data else []

        # Compile path
        workspace_path = self.cmd_pointer.workspace_path(self.cmd_pointer.settings["workspace"])

        print(indices)
        return "ok"


# Sorter function.
def _sort_mol(mol, sort, category):
    value = (mol.get(category) or {}).get(sort)
    return (value is None, value.lower() if isinstance(value, str) else value)
