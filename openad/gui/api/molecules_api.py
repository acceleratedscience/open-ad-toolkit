import json
from rdkit import Chem
from flask import Response, request

# from openad.molecules.mol_api import get_molecule_data
from openad.molecules.mol_commands import retrieve_mol
from openad.molecules.mol_functions import organize_properties
from openad.molecules.mol_functions import mol2svg, mol2sdf


class MoleculesApi:
    def __init__(self, cmd_pointer):
        self.cmd_pointer = cmd_pointer

    # Get molecule data, plus SDF and SVG.
    # Used when requesting a molecule by its identifier.
    def get_mol_data(self):
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
        inchi = data["inchi"]

        mol_rdkit = Chem.MolFromInchi(inchi)
        if mol_rdkit:
            mol_svg = mol2svg(mol_rdkit)
            mol_sdf = mol2sdf(mol_rdkit)
        else:
            mol_svg, mol_sdf = None, None

        return {"sdf": mol_sdf, "svg": mol_svg}
