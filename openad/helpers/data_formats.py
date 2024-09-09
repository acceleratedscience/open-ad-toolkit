import copy

# Template of the meta dictionary where we
# store any additional data added by user.
# This pertains future functionality for the GUI.
_META_DICT = {
    "notes": "",
    "labels:": [],
}

# Template for a small molecule dictionary.
OPENAD_MOL_DICT = {
    "name": None,
    "synonyms": [],
    "properties": {},
    "property_sources": {},
    "sources": {},
    "commments": {},
    "analysis": [],
    "enriched": False,
    # Additional data added by user.
    "meta": copy.deepcopy(_META_DICT),
}

# Template for a macromolecule dictionary.
OPENAD_MMOL_DICT = {
    # The molType determines what kind of data
    # is expected in the data dictionary.
    # Currently only "protein" is supported, but
    # in the future we may add "dna", "rna", etc.
    "molType": "",
    # Key value pairs with whatever data is available.
    "data": {},
    # Usually CIF data that can be fed to the Miew 3D viewer,
    # however this can support any other format supported by Miew.
    "data3D": "",
    # This tells Miew how to interpret the data3D field.
    # Usually this would be cif, but can also be pdb, xyz, etc
    # Please refer to the Format3D type in the openad-gui repository:
    # https://github.com/acceleratedscience/openad-gui/blob/main/src/types/index.ts
    "data3DFormat": "",
    # Additional data added by user.
    "meta": copy.deepcopy(_META_DICT),  # Additional data added by user
}


# Template of a protein
OPENAD_PROTEIN_DICT = {
    "mol_type": "protein",
    "attributes": {
        # Identification
        "idcode": "",
        "name": "",
        # Publication
        "head": "",
        "author": "",
        "release_date": "",
        "deposition_date": "",
        "keywords": "",
        "journal": "",
        "journal_reference": "",
        # Context
        "resolution": 0,
        "source": {},
        "structure_method": "",
        "structure_reference": [],
        "has_missing_residues": False,
        "missing_residues": [],
        "biomoltrans": [],
        "compound": {},
    },
    # PDB/CIF file content as string.
    # Used to generate the 3D view.
    "pdb_data": "",
    "cif_data": "",
}
