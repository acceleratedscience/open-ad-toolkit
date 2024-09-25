import copy
import gemmi
from openad.helpers.data_formats import OPENAD_MMOL_DICT
from openad.mmols.mmol_functions import parse_cif_block


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


# def cif_path2mmol(cif_path):
#     """
#     Takes the content of a .cif file and returns a macromolecule dictionary.
#     Used for opening a CIF file in the GUI.
#     """

#     # Parse the CIF file
#     cif_doc = gemmi.cif.read_file(cif_path)
#     cif_block = cif_doc.sole_block()
#     data = parse_cif_block(cif_block)

#     # Read the CIF file content
#     with open(cif_path, "r", encoding="utf-8") as f:
#         cif_data = f.read()

#     # Create the moll object
#     mmol_dict = copy.deepcopy(OPENAD_MMOL_DICT)
#     mmol_dict["molType"] = "mmol"
#     mmol_dict["data"] = data
#     mmol_dict["data3D"] = cif_data
#     mmol_dict["data3DFormat"] = "cif"

#     # Return the moll object
#     return mmol_dict, None


# def pdb_path2mmol(pdb_path):
#     """
#     Takes the content of a .pdb file and returns a macromolecule dictionary.
#     Used for opening a PDB file in the GUI.
#     """

#     # Parse the PDB file
#     parser = PDBParser(QUIET=True)
#     structure = parser.get_structure("molecule", pdb_path)
#     data = structure.header

#     # Read the PDB file content
#     with open(pdb_path, "r", encoding="utf-8") as f:
#         sdf_data = f.read()

#     # Create the moll object
#     mmol_dict = copy.deepcopy(OPENAD_MMOL_DICT)
#     mmol_dict["molType"] = "mmol"
#     mmol_dict["data"] = data
#     mmol_dict["data3D"] = sdf_data
#     mmol_dict["data3DFormat"] = "pdb"

#     # _print_all_available_pdb_data(structure, parser)

#     # Return the moll object
#     return mmol_dict, None


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
        print("cif2pdb() - No cif_data or cif_path provided")
        return None

    # Write the PDB to disk
    if dest_path:
        structure.write_pdb(dest_path)

    # Return the PDB data as a string
    else:
        return structure.make_pdb_string()


def cif2mmol(cif_data=None, cif_path=None):
    """
    Convert CIF data or file to a Moll object.

    Used for:
    - Opening a CIF file in the GUI (cif_path).
    - Opening downloaded CIF data in the GUI (cif_data).
    """

    # Parse the CIF string
    if cif_data:
        cif_doc = gemmi.cif.read_string(cif_data)

    # Load the CIF file
    elif cif_path:
        cif_doc = gemmi.cif.read_file(cif_path)
        # Read the CIF file content
        with open(cif_path, "r", encoding="utf-8") as f:
            cif_data = f.read()

    # Error handling
    else:
        print("cif2moll() - No cif_data or cif_path provided")
        return None

    # Parse the CIF data
    cif_block = cif_doc.sole_block()
    data = parse_cif_block(cif_block)

    # Create the moll object
    mmol_dict = copy.deepcopy(OPENAD_MMOL_DICT)
    mmol_dict["molType"] = "mmol"
    mmol_dict["data"] = data
    mmol_dict["data3D"] = cif_data
    mmol_dict["data3DFormat"] = "cif"

    # Return the moll object
    return mmol_dict


def pdb2mmol(pdb_data=None, pdb_path=None):
    """
    Convert PDB data or file to a Moll object.

    Used for:
    - Opening a PDB file in the GUI (pdb_path).
    - Opening downloaded PDB data in the GUI (pdb_data) **

    ** Not used because we use the CIF format when downloading
    """

    # Parse the PDB string
    if pdb_data:
        struct = gemmi.read_pdb_string(pdb_data)

    # Load the PDB file
    elif pdb_path:
        struct = gemmi.read_pdb(pdb_path)
        # Read the PDB file content
        with open(pdb_path, "r", encoding="utf-8") as f:
            pdb_data = f.read()

    # Error handling
    else:
        print("pdb2moll() - No pdb_data or pdb_path provided")
        return None

    # Parse the PDB data as CIF
    cif_doc = struct.make_mmcif_document()
    cif_block = cif_doc.sole_block()
    data = parse_cif_block(cif_block)

    # Create the moll object
    mmol_dict = copy.deepcopy(OPENAD_MMOL_DICT)
    mmol_dict["molType"] = "mmoll"
    mmol_dict["data"] = data
    mmol_dict["data3D"] = pdb_data
    mmol_dict["data3DFormat"] = "pdb"

    # Return the moll object
    return mmol_dict


#
#
#
#


# Development only
# Reveal additional PDB data we can extract, maybe for later use.
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
