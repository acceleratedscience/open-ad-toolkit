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
    # The main name of this molecule.
    "name": None,
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
