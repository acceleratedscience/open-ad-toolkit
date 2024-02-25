import json
from rdkit import Chem
from flask import request

# from openad.molecules.mol_api import get_molecule_data
from openad.molecules.mol_commands import retrieve_mol
from openad.molecules.mol_functions import organize_properties
from openad.molecules.mol_functions import mol2svg, mol2sdf


class MoleculesApi:
    def __init__(self, cmd_pointer):
        self.cmd_pointer = cmd_pointer

    def get_mol_data(self):
        data = json.loads(request.data) if request.data else {}
        identifier = data["identifier"]

        if not identifier:
            return "No identifier provided", 500

        mol = retrieve_mol(identifier)

        if mol:
            mol = organize_properties(mol)

            # Render SVG and SDF
            mol_rdkit = Chem.MolFromInchi(mol["identifiers"]["inchi"])
            if mol_rdkit:
                mol_svg = mol2svg(mol_rdkit)
                mol_sdf = mol2sdf(mol_rdkit)
            else:
                mol_svg, mol_sdf = None, None
            return {"mol": mol, "sdf": mol_sdf, "svg": mol_svg}, 200

        else:
            return "Failed to retrieve molecule", 500
