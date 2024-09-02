import copy
from Bio.PDB import PDBParser

from openad.macromolecules.mmol_functions import OPENAD_PROTEIN_DICT


def pdb_path2prot(pdb_path):
    """
    Takes the content of a .pdb file and returns a protein dictionary.
    """

    mmol_dict = copy.deepcopy(OPENAD_PROTEIN_DICT)

    # Parse the PDB file
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("molecule", pdb_path)

    # Store header data
    mmol_dict["header"] = structure.header

    # Attach the PDB file content
    with open(pdb_path, "r") as pdb_file:
        mmol_dict["pdb"] = pdb_file.read()

    print(888, mmol_dict)

    return mmol_dict, None


# def _pdb_model_to_dict(model):
#     model_dict = {
#         "model_id": model.id,
#         "chains": []
#     }
#     for chain in model:
#         chain_dict = {
#             "chain_id": chain.id,
#             "residues": []
#         }
#         for residue in chain:
#             residue_dict = {
#                 "residue_name": residue.resname,
#                 "residue_id": residue.id[1],
#                 "atoms": []
#             }
#             for atom in residue:
#                 atom_dict = {
#                     "atom_name": atom.name,
#                     "coordinates": atom.coord.tolist(),  # Convert numpy array to list
#                     "element": atom.element
#                 }
#                 residue_dict["atoms"].append(atom_dict)
#             chain_dict["residues"].append(residue_dict)
#         model_dict["chains"].append(chain_dict)
#     return model_dict
