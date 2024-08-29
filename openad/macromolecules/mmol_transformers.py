import copy
from Bio.PDB import PDBParser

from openad.macromolecules.mmol_functions import OPENAD_MMOL_DICT


def pdb_path2mol(pdb_path):
    """
    Takes the content of a .pdb file and returns a macromol dictionary.
    """

    mmol_dict = copy.deepcopy(OPENAD_MMOL_DICT)

    # Parse the PDB file
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("molecule", pdb_path)

    # Store header data
    mmol_dict.header = structure.header
    mmol_dict.models = structure.models
