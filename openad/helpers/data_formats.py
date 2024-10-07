import copy

# Template of the meta dictionary where we
# store any additional data added by user.
# This pertains future functionality for the GUI.
_META_DICT = {
    "notes": "",
    "labels:": [],
}

# Template for a small molecule dictionary.
OPENAD_SMOL_DICT = {
    # The main name of this molecule.
    "identifiers": {
        "name": None,
        "inchi": None,
        "inchikey": None,
        "canonical_smiles": None,
        "isomeric_smiles": None,
        "smiles": None,
        "molecular_formula": None,
        "cid": None,
    },
    # Alternative names this molecule is know by.
    "synonyms": [],
    # Available molecular properties.
    "properties": {
        # "molecular_weight": 0.1234
    },
    # The sources of the properties.
    # Eg. pubchem, RDKit, etc.
    "property_sources": {
        # "molecular_weight": { source: "pubchem" }
    },
    # List of analysis records for the molecule.
    # See ANALYSIS_RECORD below.
    "analysis": [],
    # Flag to indicate if the molecule data was
    # enriched with data from PubChem.
    "enriched": False,
    # Additional data added by user.
    "meta": copy.deepcopy(_META_DICT),
    #
    #
    "sources": {},  # This should go, only used once by _get_mol to update property_sources
}

# Template for an analysis record dictionary, which is
# stored under "analysis" in the small molecule dictionary.
ANALYSIS_RECORD = {
    "toolkit": None,
    "function": None,
    "smiles": None,
    "parameters": {},
    "results": [],
}

# Template for a macromolecule dictionary.
OPENAD_MMOL_DICT = {
    # The molType determines what kind of data is expected in
    # the data dictionary. Currently this will always be "mmoll",
    # but in the future we may differentiate between protein, dna, rna, etc.
    "molType": "",
    # Key value pairs with whatever data is available.
    "data": {},
    # Usually CIF or PDB data that can be fed to the Miew 3D viewer.
    "data3D": "",
    # This tells Miew how to interpret the data3D field.
    # Other than "cif" or "pdb", this can also be xyz, etc.
    # Please refer to the Format3D type in the openad-gui repository:
    # https://github.com/acceleratedscience/openad-gui/blob/main/src/types/index.ts
    "data3DFormat": "",
    # Additional data added by user.
    "meta": copy.deepcopy(_META_DICT),  # Additional data added by user
}
