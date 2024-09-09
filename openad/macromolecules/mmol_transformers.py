import copy
import gemmi
from Bio.PDB import PDBParser, MMCIFParser

from openad.macromolecules.mmol_functions import OPENAD_PROTEIN_DICT, OPENAD_MMOL_DICT


def mmol2cif(mmol_dict, path=None):
    """
    Convert a macromolecule dictionary to CIF format.
    Used to store a macromolecule as a CIF file.
    """

    # Error handling
    if not mmol_dict:
        print("mmol2pdb() - No mmol_dict provided")
        return None

    # Load the PDB data
    if mmol_dict["data3DFormat"] == "cif":
        cif_data = mmol_dict["data3D"]
        if path:
            with open(path, "w", encoding="utf-8") as f:
                f.write(cif_data)
    elif mmol_dict["data3DFormat"] == "pdb":
        cif_data = pdb2cif(mmol_dict["data3D"], dest_path=path)  # Will write to disk if path is set

    # Return the PDB data as a string
    return cif_data


def mmol2pdb(mmol_dict, path=None):
    """
    Convert a macromolecule dictionary to PDB format.
    Used to store a macromolecule as a PDB file.
    """

    # Error handling
    if not mmol_dict:
        print("mmol2pdb() - No mmol_dict provided")
        return None

    # Load the PDB data
    if mmol_dict["data3DFormat"] == "pdb":
        pdb_data = mmol_dict["data3D"]
        if path:
            with open(path, "w", encoding="utf-8") as f:
                f.write(pdb_data)

    elif mmol_dict["data3DFormat"] == "cif":
        pdb_data = cif2pdb(mmol_dict["data3D"], dest_path=path)  # Will write to disk if path is set

    # Return the PDB data as a string
    return pdb_data


def cif_path2mmol(cif_path):
    """
    Takes the content of a .cif file and returns a macromolecule dictionary.
    Used for opening a CIF file in the GUI.
    """

    mmol_dict = copy.deepcopy(OPENAD_PROTEIN_DICT)

    # Parse the PDB file
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("molecule", cif_path)
    data = structure.header

    # Load the CIF file content
    with open(cif_path, "r", encoding="utf-8") as f:
        cif_data = f.read()

    # Fill in the data
    mmol_dict["molType"] = "protein"
    mmol_dict["data"] = data
    mmol_dict["data3D"] = cif_data
    mmol_dict["data3DFormat"] = "cif"

    return mmol_dict, None


def pdb_path2mmol(pdb_path):
    """
    Takes the content of a .pdb file and returns a macromolecule dictionary.
    Used for opening a PDB file in the GUI.
    """

    mmol_dict = copy.deepcopy(OPENAD_MMOL_DICT)

    # Parse the PDB file
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("molecule", pdb_path)
    data = structure.header

    # Read the PDB file content
    with open(pdb_path, "r", encoding="utf-8") as f:
        sdf_data = f.read()

    # Fill in the data
    mmol_dict["molType"] = "protein"
    mmol_dict["data"] = data
    mmol_dict["data3D"] = sdf_data
    mmol_dict["data3DFormat"] = "pdb"

    # _print_all_available_pdb_data(structure, parser)

    return mmol_dict, None


def pdb2cif(pdb_data=None, pdb_path=None, dest_path=None):
    """
    Convert PDB path/data to CIF format.
    """

    # Note: converting a PDF into CIF format does not seem to work well.
    # Neither gemmi or biopython manages to preserve the data in a way
    # that it stays compatible with our Miew viewer. Because of this,
    # we've disable converting from PDB to CIF on the GUI side, but we
    # keep this function in place in case this changes in the future.

    print("pdb2cif() is disabled.")
    return None

    # Parse the PDB string
    if pdb_data:
        structure = gemmi.read_pdb_string(pdb_data)

    # Load the PDB file
    elif pdb_path:
        structure = gemmi.read_structure(pdb_path)

    # Error handling
    else:
        print("pdb2cif() - No pdb_path or pdb_data provided")
        return None

    # See https://gemmi.readthedocs.io/en/latest/mol.html#entity
    structure.setup_entities()
    structure.assign_label_seq_id()

    # Write the CIF to disk
    if dest_path:
        structure.make_mmcif_document().write_file(dest_path)

    # Return the CIF data as a string
    else:
        return structure.make_mmcif_block()


def cif2pdb(cif_data=None, cif_path=None, dest_path=None):
    """
    Convert CIF path to PDB format.
    """

    # Parse the CIF string
    if cif_data:
        cif_doc = gemmi.cif.read_string(cif_data)
        cif_block = cif_doc.sole_block()
        structure = gemmi.make_structure_from_block(cif_block)

    # Load the CIF file
    elif cif_path:
        structure = gemmi.read_structure(cif_path)

    # Error handling
    else:
        print("cif2pdb() - No cif_path or cif_data provided")
        return None

    # Write the PDB to disk
    if dest_path:
        # structure = gemmi.make_structure(cif_doc)
        # cif_doc.write_file(dest_path)

        structure.write_pdb(dest_path)
        # cif_doc = structure.make_mmcif_document()
        # cif_doc.write_file(dest_path)

    # Return the PDB data as a string
    else:
        return structure.make_pdb_string()


#
#
#
#


# Development only
# Reveal additional PDB data to extract - currently not used.
def _print_all_available_pdb_data(structure, parser):
    struct_id = structure.id
    level = structure.level
    child_dict = structure.child_dict
    trailer = parser.get_trailer()
    models = structure.get_models()
    chains = structure.get_chains()
    residues = structure.get_residues()
    atoms = structure.get_atoms()

    # Print id, level
    print("\n\nId:\n", struct_id)
    print("\n\nLevel:\n", level)

    # Print child_dict
    print("\n\nChild dict:")
    for key, value in child_dict.items():
        print(f"Key: {key}, Value: {value}")

    # Print trailer
    print("\n\nTrailer:\n", trailer)

    # Print models
    print("\n\nModels:")
    for model in models:
        print(f"Model ID: {model.id}")
        for chain in model:
            print(f"  Chain ID: {chain.id}")
            for residue in chain:
                print(f"    Residue: {residue.resname} {residue.id}")
                for atom in residue:
                    print(f"      Atom: {atom.name} {atom.coord}")

    # Print chains
    print("\n\nChains:")
    for chain in chains:
        print(f"Chain ID: {chain.id}")
        for residue in chain:
            print(f"  Residue: {residue.resname} {residue.id}")
            for atom in residue:
                print(f"    Atom: {atom.name} {atom.coord}")

    # Print residues
    print("\n\nResidues:")
    for residue in residues:
        print(f"Residue: {residue.resname} {residue.id}")
        for atom in residue:
            print(f"  Atom: {atom.name} {atom.coord}")

    # Print atoms
    print("\n\nAtoms:")
    for atom in atoms:
        print(f"Atom: {atom.name} {atom.coord}")
