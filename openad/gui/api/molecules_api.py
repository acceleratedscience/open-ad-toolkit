import json
import math
from rdkit import Chem
from flask import Response, request

# from openad.molecules.mol_api import get_molecule_data
from openad.molecules.mol_commands import retrieve_mol
from openad.molecules.mol_functions import organize_properties
from openad.molecules.mol_functions import mol2svg, mol2sdf

from openad.helpers.json_decimal_encoder import DecimalEncoder
from openad.helpers.files import open_file


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
        path = data["path"] if "path" in data else ""
        query = data["query"] if "query" in data else ""
        page = data["page"] if "page" in data else 1
        page_size = data["pageSize"] if "pageSize" in data else 48

        # Compile path
        workspace_path = self.cmd_pointer.workspace_path(self.cmd_pointer.settings["workspace"])
        file_path = workspace_path + "/" + path

        # Read file
        molset, err_code = open_file(file_path, return_err="code", as_string=False)

        # Return error
        if err_code:
            return None

        # Search all identifiers and properties for the query.
        results = []
        for mol in molset:
            found = False
            for key in mol["identifiers"]:
                if query in str(mol["identifiers"][key]):
                    results.append(mol)
                    break
            if found:
                continue
            for key in mol["properties"]:
                if query in str(mol["properties"][key]):
                    results.append(mol)
                    break

        total = len(results)
        page_total = math.ceil(total / page_size)

        skip = 48 * (page - 1)
        result_page = json.dumps(results[skip : skip + 48], cls=DecimalEncoder)

        data = {
            "total": len(results),
            "mols": result_page,
            "page": page,
            "pageSize": page_size,
            "pageTotal": page_total,
        }

        return data
