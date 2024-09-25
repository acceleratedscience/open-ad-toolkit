"""
Molecule management functions
"""

import os
import ast
import time
import json
import shutil
import pandas
import asyncio
import aiofiles
import pubchempy as pcy
from copy import deepcopy
from datetime import datetime
from rdkit import Chem, rdBase
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.Descriptors import MolWt, ExactMolWt

# OpenAD imports
from openad.helpers.output_msgs import msg
from openad.helpers.output import output_error, output_warning, output_success
from openad.helpers.files import open_file
from openad.helpers.general import confirm_prompt, pretty_date
from openad.helpers.spinner import spinner
from openad.helpers.json_decimal_encoder import DecimalEncoder
from openad.helpers.data_formats import OPENAD_SMOL_DICT

# Suppress RDKit logs
rdBase.BlockLogs()  # pylint: disable=c-extension-no-member

# PubChem accepted molecule identifiers
PCY_IDFR = {
    "name": "name",
    "smiles": "smiles",
    "inchi": "inchi",
    "inchikey": "inchikey",
    "cid": "cid",
    "formula": "formula",
}
MOL_PROPERTY_SOURCES = {
    "Log P-XLogP3-AA": "xlogp",
    "Log P-XLogP3": "xlogp",
    "SMILES-Isomeric": "isomeric_smiles",
    "SMILES-Canonical": "canonical_smiles",
    "Molecular Weight": "molecular_weight",
    "Compound Complexity": "complexity",
    "Count-Rotatable Bond": "rotatable_bond_count",
    "Compound-Canonicalized": "complexity",
    "Count-Hydrogen Bond Acceptor": "h_bond_acceptor_count",
    "Count-Hydrogen Bond Donor": "h_bond_donor_count",
    "IUPAC Name-Preferred": "iupac_name",
    "Fingerprint-SubStructure Keys": "",
    "InChI-Standard": "inchi",
    "InChIKey-Standard": "inchikey",
    "Mass-Exact": "exact_mass",
    "Weight-MonoIsotopic": "monoisotopic_mass",
    "Molecular Formula": "molecular_formula",
    "Topological-Polar Surface Area": "tpsa",
}

# Molecule properties with machine data which we omit
# when displaying the molecule properties in the CLI
SMOL_PROPS_IGNORE = [
    "atoms",
    "bonds",
    "cactvs_fingerprint",
    "elements",
    "fingerprint",
    "record",
]

# @@Todo: retire thsi in favor of SMOL_PROPS_IGNORE?
SMOL_PROPERTIES = sorted(
    [
        "atom_stereo_count",
        "bond_stereo_count",
        "canonical_smiles",
        "charge",
        "cid",
        "complexity",
        "conformer_count_3d",
        "conformer_id_3d",
        "conformer_model_rmsd_3d",
        "conformer_rmsd_3d",
        "coordinate_type",
        "covalent_unit_count",
        "defined_atom_stereo_count",
        "defined_bond_stereo_count",
        "effective_rotor_count_3d",
        "exact_mass",
        "feature_acceptor_count_3d",
        "feature_anion_count_3d",
        "feature_cation_count_3d",
        "feature_count_3d",
        "feature_donor_count_3d",
        "feature_hydrophobe_count_3d",
        "feature_ring_count_3d",
        "h_bond_acceptor_count",
        "h_bond_donor_count",
        "heavy_atom_count",
        "inchi",
        "inchikey",
        "isomeric_smiles",
        "isotope_atom_count",
        "iupac_name",
        "mmff94_energy_3d",
        "mmff94_partial_charges_3d",
        "molecular_formula",
        "molecular_weight_exact",
        "molecular_weight",
        "monoisotopic_mass",
        "multipoles_3d",
        "multipoles_3d",
        "pharmacophore_features_3d",
        "pharmacophore_features_3d",
        "rotatable_bond_count",
        "sol_classification",
        "sol",
        "tpsa",
        "undefined_atom_stereo_count",
        "undefined_bond_stereo_count",
        "volume_3d",
        "x_steric_quadrupole_3d",
        "xlogp",
        "y_steric_quadrupole_3d",
        "z_steric_quadrupole_3d",
    ]
)
INPUT_MAPPINGS = {
    "NAME": "chemical_name",
    "xlogp3": "xlogp",
    "molecular weight": "molecular_weight",
    "complexity": "complexity",
    "rotatable bond count": "rotatable_bond_count",
    "hydrogen bond acceptor count": "h_bond_acceptor_count",
    "hydrogen bond donor count": "h_bond_donor_count",
    "exact mass": "exact_mass",
    "monoisotopic mass": "monoisotopic_mass",
    "topological polar surface area": "tpsa",
    "heavy atom count": "heavy_atom_count",
    "formal charge": "formal_charge",
    "isotope atom count": "isotope_atom_count",
    "defined atom stereocenter count": "defined_atom_stereo_count",
    "undefined atom stereocenter count": "undefined_atom_stereo_count",
    "covalently-bonded unit count": "covalent_unit_count",
    "compound is canonicalized": "compound_canonicalized",
    "SOL_classification": "sol_classification",
    "SOL": "sol",
}

mol_name_cache = {}

############################################################
# region - Lookup / creation


def find_smol(cmd_pointer: object, identifier: str, name: str = None, basic: bool = False) -> dict | None:
    """
    Find a molecule across the available resources.
    First we check our working set, then PubChem, then RDKit.

    Parameters
    ----------
    cmd_pointer: object
        The command pointer object.
    identifier: str
        The molecule identifier to search for.
        Valid inputs: InChI, SMILES, InChIKey, name, CID.
    name: str
        Optional name for the molecule.
    basic: bool
        If True, create a basic molecule dict with RDKit (no API calls).

    Returns
    -------
    smol: dict
        The OpenAD small molecule dictionary if found, otherwise None.
    """

    # Look for molecule in the working set
    smol = get_smol_from_mws(cmd_pointer, identifier)

    # Look for molecule in memory
    if not smol:
        smol = get_smol_from_memory(cmd_pointer, identifier)

    # Look for molecule on PubChem
    if not smol and not basic:
        smol = get_smol_from_pubchem(identifier)

    # Try creating molecule object with RDKit.
    if not smol:
        smol = new_smol(identifier, name=name)

    # Fail - invalid.
    if not smol:
        output_error(["Unable to identify molecule", f"Identifier: {identifier}"], return_val=False)
        if basic is True:
            output_error(["Invalid InChI or SMILES string", f"Identifier: {identifier}"], return_val=False)

    return smol


def get_smol_from_mws(cmd_pointer: object, identifier: str, ignore_synonyms: bool = False) -> dict | None:
    """
    Retrieve a molecule from the molecule working set.

    Parameters
    ----------
    cmd_pointer: object
        The command pointer object.
    identifier: str
        The molecule identifier to search for.
        Valid inputs: InChI, SMILES, InChIKey, name, CID.
    ignore_synonyms: bool
        If True, ignore synonyms in the search.

    Returns
    -------
    dict
        The OpenAD smol dictionary if found, otherwise None.
    """

    smol = get_smol_from_list(identifier, cmd_pointer.molecule_list, ignore_synonyms=ignore_synonyms)
    if smol is not None:
        return deepcopy(smol)
    return None


def get_smol_from_memory(cmd_pointer: object, identifier: str) -> dict | None:
    """
    Retrieve a molecule from memory.

    Parameters
    ----------
    cmd_pointer: object
        The command pointer object.
    identifier: str
        The molecule identifier to search for.
        Valid inputs: InChI, SMILES, InChIKey, name, CID.
    """

    if cmd_pointer.last_external_molecule is not None:
        return get_smol_from_list(identifier, [cmd_pointer.last_external_molecule])
    return None


def get_smol_from_pubchem(identifier: str) -> dict | None:
    """
    Fetch small molecule from PubChem.

    Parameters
    ----------
    identifier: str
        The small molecule identifier to search for.
        Valid inputs: InChI, SMILES, InChIKey, name, CID.
    """

    smol = _get_pubchem_compound(identifier, PCY_IDFR["inchi"])
    if not smol:
        smol = _get_pubchem_compound(identifier, PCY_IDFR["smiles"])
    if not smol:
        smol = _get_pubchem_compound(identifier, PCY_IDFR["inchikey"])
    if not smol:
        smol = _get_pubchem_compound(identifier, PCY_IDFR["name"])
    if not smol:
        smol = _get_pubchem_compound(identifier, PCY_IDFR["cid"])
    # # Commented out until getting no time outs from pubchem.
    # if not smol:
    #     smol = _get_pubchem_compound(identifier, PCY_IDFR["formula"])

    if smol:
        return _sep_identifiers_from_properties(smol)

    return None


# region--local
def _get_pubchem_compound(identifier: str, identifier_type: str) -> dict | None:
    """
    Fetch small molecule from PubChem based on an identifier.

    Parameters
    ----------
    identifier: str
        The small molecule identifier to search for.
        Valid inputs: InChI, SMILES, InChIKey, name, CID.
    identifier_type: str
        The type of identifier to search for (see PCY_IDFR).
    """

    mol_pcy = None

    try:
        # Find molecule on PubChem
        compounds = pcy.get_compounds(identifier, identifier_type)
        if len(compounds) == 0:
            return None
        else:
            mol_pcy = compounds[0].to_dict()

        # Create OpenAD smol dict
        if mol_pcy:
            smol = deepcopy(OPENAD_SMOL_DICT)
            smol = _add_pcy_data(smol, mol_pcy, identifier, identifier_type)
            return smol

    except Exception:  # pylint: disable=broad-exception-caught
        pass

        # # Keep here for debugging
        # output_error(
        #     [
        #         "Error _get_pubchem_compound()",
        #         f"identifier: {identifier}",
        #         f"identifier_type: {identifier_type}",
        #     ]
        # )

    return None


def _add_pcy_data(smol, smol_pcy, identifier, identifier_type):
    """
    Add PubChem molecule data to the OpenAD molecule dict.

    Parameters
    ----------
    smol: dict
        The OpenAD small molecule dictionary to which we add the PubChem data.
    mol_pcy: dict
        The PubChem molecule data.
    identifier: str
        The molecule identifier.
    identifier_type: str
        The type of identifier to search for (see PCY_IDFR).
    """

    smol["enriched"] = True

    # Add synonyms
    synonyms = pcy.get_synonyms(smol_pcy["iupac_name"], "name")
    smol["synonyms"] = synonyms[0].get("Synonym") if synonyms and len(synonyms) > 0 else []

    # Add name
    if identifier_type == PCY_IDFR["name"]:
        smol["name"] = identifier
    elif len(smol.get("synonyms", [])) > 1:
        smol["name"] = smol["synonyms"][0]
    else:
        smol["name"] = smol_pcy["iupac_name"]

    # Add properties
    smol["properties"] = smol_pcy

    # Add canonical smiles
    if identifier_type == PCY_IDFR["smiles"]:
        # fmt: off
        smol["properties"]["canonical_smiles"] = canonicalize(identifier)  # pylint: disable=no-member
        # fmt: on

    # Loop through PubChem properties and update our own property_sources
    # when PubChem has its own 3rd party source for a property.
    # - - -
    # For example the property_sources["iupac_name"] value when looking up dopamine:
    # - Before: {"source": "pubchem"}
    # - After: { 'label': 'IUPAC Name', 'name': 'Preferred', 'datatype': 1, 'version': '2.7.0',
    #            'software': 'Lexichem TK', 'source': 'OpenEye Scientific Software', 'release': '2021.10.14'}
    for x in SMOL_PROPERTIES:
        smol["property_sources"][x] = {"source": "PubChem"}
        for prop_name, prop_name_key in MOL_PROPERTY_SOURCES.items():
            if prop_name_key == x:
                if len(prop_name.split("-")) > 0:
                    for y in smol_pcy["record"]["props"]:
                        if "label" not in y["urn"]:
                            pass
                        elif y["urn"]["label"] == prop_name.split("-", maxsplit=1)[0] and "name" not in y["urn"]:
                            smol["property_sources"][x] = y["urn"]
                        elif (
                            y["urn"]["label"] == prop_name.split("-", maxsplit=1)[0]
                            and y["urn"]["name"] == prop_name.split("-", maxsplit=2)[1]
                        ):
                            smol["property_sources"][x] = y["urn"]

    return smol


def _sep_identifiers_from_properties(smol: dict) -> dict:
    """
    Separate molecules identifiers from properties.

    This is the final step when processing external molecule data
    from a file like MDL, SDF, CSV, etc. or from an API call, so
    a molecule is in the correct OpenAD format.

    Parameters:
    -----------
    smol: dict
        The molecule object to modify.
    """

    # Move all identifiers to the identifiers key.
    smol["identifiers"] = _get_identifiers(smol)

    # Remove identifiers from properties.
    # Create a lowercase version of the properties dictionary
    # so we can scan for properties in a case-insensitive way.
    molIdfrs = {k.lower(): v for k, v in smol["identifiers"].items()}
    for prop in list(smol["properties"]):
        if prop.lower() in molIdfrs:
            del smol["properties"][prop]

    return smol


def _get_identifiers(smol: dict) -> dict:
    """
    Pull the identifiers from a molecule.
    """

    # In case this smol has the identifiers already separated.
    if smol.get("identifiers"):
        return smol["identifiers"]

    identifier_dict = {}

    # Create a lowercase version of the properties dictionary
    # so we can scan for properties in a case-insensitive way.
    smolProps = {k.lower(): v for k, v in smol["properties"].items()}

    identifier_dict["name"] = smolProps.get("name")
    identifier_dict["inchi"] = smolProps.get("inchi")
    identifier_dict["inchikey"] = smolProps.get("inchikey")
    identifier_dict["canonical_smiles"] = smolProps.get("canonical_smiles")
    identifier_dict["isomeric_smiles"] = smolProps.get("isomeric_smiles")
    identifier_dict["smiles"] = smolProps.get("smiles")
    identifier_dict["molecular_formula"] = smolProps.get("molecular_formula")
    identifier_dict["cid"] = smolProps.get("cid")

    return identifier_dict


# endregion


# Takes any identifier and creates a minimal molecule object
# on the fly using RDKit, without relying on PubChem or API calls.
# - - -
# Note: in our molecule data structure, identifiers are all stored
# under properties. The GUI and possibly other parts of the
# application consume a slightly modified format, where identifiers
# are stored separately from properties. This is a cleaner / more
# correct way of organizing the molecule object, since identifiers
# are not properties, and they are treated differently (eg. no sort).
# But we can't change the main molecule datastructure without creating
# a formatter to ensure backwards compatibilty, so for now you can
# use molformat_v2() to convert the molecule object to the new format.
# - - -
# It is recommended to start using the new format elsewhere in new code,
# so we'll have less to refactor once we switch to the new format.
def new_smol(inchi_or_smiles: str = None, mol_rdkit: Mol = None, name: str = None):
    """
    Create a basic molecule object without relying on API calls.

    Parameters
    ----------
    inchi_or_smiles: str
        Source option A: An InChI or SMILES identifier.
    mol_rdkit: RDKit ROMol
        Source option B: An RDKit molecule object.
    name: str
        Optional name for the molecule.
    """

    smol = deepcopy(OPENAD_SMOL_DICT)
    timestamp = pretty_date()
    prop_src = {"source": "RDKit", "date": timestamp}

    # Create RDKit molecule object
    if mol_rdkit is None:
        try:
            mol_rdkit = Chem.MolFromInchi(inchi_or_smiles)
            if not mol_rdkit:
                mol_rdkit = Chem.MolFromSmiles(inchi_or_smiles)  # pylint: disable=no-member
            if not mol_rdkit:
                mol_rdkit = Chem.MolFromInchi("InChI=1S/" + inchi_or_smiles)
            if not mol_rdkit:
                return None
        except Exception:  # pylint: disable=broad-exception-caught
            return None

    # Parse properties
    props = mol_rdkit.GetPropsAsDict()
    for key in props:
        smol["properties"][key] = props[key]
        smol["property_sources"][key] = prop_src

    # Figure out the best name
    if name is None:
        name = smol["properties"].get("name", None)

    # fmt: off
    # Store identifiers
    smol["name"] = name
    smol["identifiers"]["inchi"] = Chem.MolToInchi(mol_rdkit)
    smol["identifiers"]["inchikey"] = Chem.inchi.InchiToInchiKey(smol["properties"]["inchi"])
    smol["identifiers"]["canonical_smiles"] = Chem.MolToSmiles(mol_rdkit)  # pylint: disable=no-member
    smol["identifiers"]["isomeric_smiles"] = Chem.MolToSmiles(mol_rdkit, isomericSmiles=True)  # pylint: disable=no-member
    smol["identifiers"]["molecular_formula"] = Chem.rdMolDescriptors.CalcMolFormula(mol_rdkit) # pylint: disable=c-extension-no-member
    smol["properties"]["molecular_weight"] = MolWt(mol_rdkit)
    
    # Store property sources
    smol["property_sources"]["name"] = prop_src
    smol["property_sources"]["inchi"] = prop_src
    smol["property_sources"]["inchikey"] = prop_src
    smol["property_sources"]["canonical_smiles"] = prop_src
    smol["property_sources"]["isomeric_smiles"] = prop_src
    smol["property_sources"]["molecular_formula"] = prop_src
    smol["property_sources"]["molecular_weight"] = prop_src

    # Disabled this for consistency, because molecular_weight_exact is not included in PubChem data.
    # See: pcy.Compound.from_cid(cid).to_dict()
    # mol["properties"]["molecular_weight_exact"] = ExactMolWt(mol_rdkit)
    # mol["property_sources"]["molecular_weight_exact"] = {"source": "RDKit", "date": timestamp}
    # fmt: on

    # So the UI can recognize when a molecule has been enriched.
    smol["enriched"] = False

    return smol


def get_human_properties(smol: dict) -> dict:
    """
    Pull subset of properties from a molecule,
    ignoring all the data meant for machines.

    Used to display a molecules's properties in the CLI.

    Parameters
    ----------
    smol: dict
        The OpenAD small molecule dictionary to extract properties from.

    Returns
    -------
    dict
        The extracted properties.
    """

    props = {}
    for prop in smol["properties"]:
        if prop not in SMOL_PROPS_IGNORE:
            props[prop] = smol["properties"][prop]
    return props


def get_molset_mols(path_absolute: str):
    """
    Return the list of molecules from a molset file,
    with an index added to each molecule.

    Parameters
    ----------
    path_absolute: str
        The absolute path to the molset file.

    Returns
    -------
    molset
        The list of molecules.
    err_code
        The open_file error code if an error occurred.
    """

    # Read file contents
    molset, err_code = open_file(path_absolute, return_err="code")

    return molset, err_code


# endregion

############################################################
# region - Validation


def valid_smiles(smiles: str) -> bool:
    """
    Verify if a string is valid SMILES definition.

    Parameters
    ----------
    smiles: str
        The SMILES string to validate
    """

    try:
        m = Chem.MolFromSmiles(smiles, sanitize=False)  # pylint: disable=no-member
    except Exception:  # pylint: disable=broad-exception-caught
        return False

    if m is None:
        return False
    else:
        try:
            Chem.SanitizeMol(m)  # pylint: disable=no-member
        except Exception:  # pylint: disable=broad-exception-caught
            return False

    return True


def valid_inchi(inchi: str) -> bool:
    """
    Verify if a string is valid InChI definition.

    Parameters
    ----------
    inchi: str
        The InChI string to validate
    """

    try:
        m = Chem.inchi.InchiToInchiKey(inchi)
    except Exception:  # pylint: disable=broad-exception-caught
        return False

    if m is None:
        return False
    else:
        return True


# endregion

############################################################
# region - Utility


def flatten_smol(smol: dict) -> dict:
    """
    Flatten identifiers, properties, name and synonyms onto
    the root of the molecule dictionary, so the molecule can
    be displayed in a table.

    Parameters
    ----------
    smol: dict
        The OpenAD small molecule dictionary to flatten.

    Returns
    -------
    dict
        The flattened molecule dictionary.
    """
    if smol is None:
        return None

    smol_flat = {**smol.get("identifiers"), **smol.get("properties")}
    smol_flat["name"] = smol.get("name", None)
    smol_flat["synonyms"] = "\n".join(smol.get("synonyms", []))
    return smol_flat


def canonicalize(smiles: str) -> str | None:
    """
    Turn any SMILES into its canonical equivalent per RDKit.

    Parameters
    ----------
    smiles: str
        The SMILES string to canonicalize
    """
    if not smiles:
        return None

    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles), isomericSmiles=True)  # pylint: disable=no-member


def df_has_molecules(df: pandas.DataFrame) -> bool:
    """
    Check if a dataframe has molecules.
    """

    cols_lowercase = [col.lower() for col in df.columns]
    if "smiles" in cols_lowercase:
        return True
    elif "inchi" in cols_lowercase:
        return True
    else:
        return False


def get_smol_from_list(identifier, molset, ignore_synonyms=False):
    """
    Scan a molset for a given identifier.
    Returns the molecule dict if found, otherwise None.

    Parameters
    ----------
    identifier: str
        The molecule identifier to search for.
        Valid inputs: InChI, SMILES, InChIKey, name, CID.
    molset: list
        A list of OpenAD small molecule objects.
    ignore_synonyms: bool
        If True, ignore synonyms in the search.
        This is only used when renaming a molecule to one of its synonyms.
        Without it, the search would return the original molecule and abort.

    """

    for openad_mol in molset:
        identifiers_dict = openad_mol.get("identifiers")
        synonyms = openad_mol.get("synonyms", [])
        if not identifiers_dict:
            output_error("get_smol_from_list(): Invalid molset input", return_val=False)
            return None

        # Name match
        if identifier.upper() == identifiers_dict.get("name", "").upper():
            return openad_mol

        # CID match
        if int(identifier) == int(identifiers_dict.get("cid", "")):
            return openad_mol

        # InChI match
        if identifier == identifiers_dict.get("inchi"):
            return openad_mol

        # InChIKey match
        if identifier == identifiers_dict.get("inchikey"):
            return openad_mol

        # SMILES match
        try:
            idfr_canonical = canonicalize(identifier)
            if idfr_canonical == canonicalize(identifiers_dict.get("canonical_smiles", None)):
                return openad_mol
            elif idfr_canonical == canonicalize(identifiers_dict.get("icomeric_smiles", None)):
                return openad_mol
            elif idfr_canonical == canonicalize(identifiers_dict.get("smiles", None)):
                return openad_mol
        except Exception:  # pylint: disable=broad-except
            pass

        # Synonym match
        if not ignore_synonyms:
            for syn in synonyms:
                if identifier.upper() == syn.upper():
                    return openad_mol

    # Fail
    return None


def get_best_available_identifier(smol: dict) -> tuple:
    """
    Get whatever identifier is available from a molecule.

    Parameters
    ----------
    smol: dict
        The OpenAD small molecule dictionary.

    Returns
    -------
    tuple
        The identifier type and the identifier string.
    """

    identifiers_dict = smol.get("identifiers", {})
    if not identifiers_dict:
        print(3333)
        return None, None

    # InChI
    inchi = identifiers_dict.get("inchi")
    if inchi:
        return "inchi", inchi

    # Isomeric SMILES
    isomeric_smiles = identifiers_dict.get("isomeric_smiles")
    if isomeric_smiles:
        return "isomeric_smiles", isomeric_smiles

    # Canonical SMILES
    canonical_smiles = identifiers_dict.get("canonical_smiles")
    if canonical_smiles:
        return "canonical_smiles", canonical_smiles

    # SMILES
    smiles = identifiers_dict.get("smiles")
    if smiles:
        return "smiles", smiles

    # InChIKey
    inchikey = identifiers_dict.get("inchikey")
    if inchikey:
        return "inchikey", inchikey

    # Name
    name = identifiers_dict.get("name")
    if name:
        return "name", name

    # CID
    cid = identifiers_dict.get("cid")
    if cid:
        return "cid", cid

    # Fail
    return None, None


def get_best_available_smiles(smol: dict) -> str | None:
    """
    Get the best available SMILES string from a molecule.

    Parameters
    ----------
    smol: dict
        The OpenAD small molecule dictionary.

    Returns
    -------
    str
        The best available SMILES string.
    """

    identifiers_dict = smol.get("identifiers")

    # Isomeric SMILES
    isomeric_smiles = identifiers_dict.get("isomeric_smiles")
    if isomeric_smiles:
        return isomeric_smiles

    # Canonical SMILES
    canonical_smiles = identifiers_dict.get("canonical_smiles")
    if canonical_smiles:
        return canonical_smiles

    # SMILES
    smiles = identifiers_dict.get("smiles")
    if smiles:
        return smiles

    # Fail
    return None


def normalize_mol_df(mol_df: pandas.DataFrame, batch: bool = False) -> pandas.DataFrame:
    """
    Normalize the column names of a molecule dataframe
    """

    has_name = False
    contains_name = None

    for col in mol_df.columns:
        # Find the name column.
        if str(col.upper()) == "NAME" or str(col.lower()) == "chemical_name":
            has_name = True
        if contains_name is None and "NAME" in str(col.upper()):
            contains_name = col
        if contains_name is None and "CHEMICAL_NAME" in str(col.upper()):
            contains_name = col

        # Normalize any columns we'll be referring to later.
        if str(col.upper()) == "SMILES":
            mol_df.rename(columns={col: "SMILES"}, inplace=True)
        if str(col.upper()) == "ROMOL":
            mol_df.rename(columns={col: "ROMol"}, inplace=True)
        if str(col.upper()) == "IMG":
            mol_df.rename(columns={col: "IMG"}, inplace=True)
        if col in INPUT_MAPPINGS:
            mol_df.rename(columns={col: INPUT_MAPPINGS[col]}, inplace=True)

    # Normalize name column.
    if has_name is False and contains_name is not None:
        mol_df.rename(columns={contains_name: "NAME"}, inplace=True)

    # Add names when missing.
    try:
        if has_name is False:
            output_warning(msg("no_m2g_name_column"))  # @@todo clean

            if not batch:
                spinner.start("Downloading names")

            mol_df["NAME"] = "unknown"
            for col in mol_df.itertuples():
                mol_df.loc[col.Index, "NAME"] = _smiles_to_iupac(mol_df.loc[col.Index, "SMILES"])

            if not batch:
                spinner.succeed("Names downloaded")
                spinner.start()
                spinner.stop()

    except Exception as err:  # pylint: disable=broad-exception-caught
        spinner.fail("There was an issue loading the molecule names.")
        spinner.start()
        spinner.stop()
        print(err)

    return mol_df


def _smiles_to_iupac(smiles):
    """
    Get the official IUPAC(*) name of a molecules based on its SMILES.

    (*) International Union of Pure and Applied Chemistry
    """

    if smiles in mol_name_cache:
        return mol_name_cache[smiles]
    try:
        compounds = pcy.get_compounds(smiles, namespace="smiles")
        match = compounds[0]
        mol_name_cache[smiles] = str(match)
    except Exception:  # pylint: disable=broad-exception-caught
        match = smiles
    return str(match)


# endregion

############################################################
# region - GUI operations


def assemble_cache_path(cmd_pointer: object, file_type: str, cache_id: str) -> str:
    """
    Compile the file path to a cached working copy of a file.

    Parameters
    ----------
    cmd_pointer: object
        The command pointer object, used to fetch the workspace path.
    file_type: 'molset'
        The type of file, used to name the cache file. For now only molset.
    cache_id: str
        The cache ID of the file.
    """

    workspace_path = cmd_pointer.workspace_path()
    return f"{workspace_path}/._openad/wc_cache/{file_type}-{cache_id}.json"


def create_molset_cache_file(cmd_pointer: object, molset: dict = None, path_absolute: str = None) -> str:
    """
    Store molset as a cached file so we can manipulate it in the GUI,
    return its cache_id so we can reference it later.

    Parameters
    ----------
    cmd_pointer: object
        The command pointer object.
    molset: dict
        The molset to cache.
    path_absolute: str
        The absolute path to the molset file, if it exists.

    Returns
    -------
    str
        The cache ID of the molset file.
    """

    cache_id = str(int(time.time() * 1000))
    cache_path = assemble_cache_path(cmd_pointer, "molset", cache_id)

    # Creaste the /._openad/wc_cache directory if it doesn't exist.
    os.makedirs(os.path.dirname(cache_path), exist_ok=True)

    # For JSON files, we can simply copy the original file (fast).
    if path_absolute:
        # timeit("copy_file")
        shutil.copy(path_absolute, cache_path)
        # timeit("copy_file", True)

        # Add indices to molecules in our working copy,
        # without blocking the thread.
        # timeit("index_wc")
        index_molset_file_async(cache_path)
        # timeit("index_wc", True)

    # For all other cases, i.e. other file formats or molset data from memory,
    # we need to write the molset object to disk (slow).
    else:
        # timeit("write_cache")
        with open(cache_path, "w", encoding="utf-8") as f:
            json.dump(molset, f)
        # timeit("write_cache", True)

    return cache_id


def read_molset_from_cache(cmd_pointer: object, cache_id: str) -> dict:
    """
    Read a cached molset file from disk.

    Parameters
    ----------
    cmd_pointer: object
        The command pointer object, used to fetch the workspace path.
    cache_id: str
        The cache ID of the molset file.

    Returns
    -------
    dict
        The molset object.
    """

    # Read file from cache.
    cache_path = assemble_cache_path(cmd_pointer, "molset", cache_id)
    molset, err_code = get_molset_mols(cache_path)

    # Return error
    if err_code:
        raise ValueError(f"Failed to read molset {cache_id} from cache")
    else:
        return molset


# endregion

############################################################
# region - Molecule working set

# Note: get_smol_from_mws() is listed under Lookup / creation.


def mws_add(cmd_pointer: object, smol: dict, force: bool = False, suppress: bool = False):
    """
    Add a molecule to your molecule working set.

    Parameters
    ----------
    cmd_pointer: object
        The command pointer object.
    smol: dict
        The OpenAD molecule object to add.
    force: bool
        If True, add without confirming.
    suppress: bool
        If True, suppress success output.
    """

    if not smol:
        output_error("No molecule provided", return_val=False)
        return False

    # Fail - already in list.
    if get_smol_from_mws(cmd_pointer, smol["properties"]["canonical_smiles"]) is not None:
        output_error("Molecule already in list: " + smol["properties"]["canonical_smiles"], return_val=False)
        return True

    # Name
    name = smol["name"]

    # Add function
    def _add_mol():
        cmd_pointer.molecule_list.append(smol.copy())
        if suppress is False:
            output_success(f"Molecule <yellow>{name}</yellow> was added", pad=0, return_val=False)
        return True

    # Add without confirming.
    if force:
        return _add_mol()

    # Confirm before adding.
    smiles = get_best_available_smiles(smol)
    smiles_str = f" <reset>{smiles}</reset>" if smiles else ""
    if confirm_prompt(f"Add molecule <green>{name}</green>{smiles_str} to your molecules working list?"):
        return _add_mol()
    else:
        output_error(f"Molecule <yellow>{name}</yellow> was not added", pad=0, return_val=False)
        return False


def mws_remove(cmd_pointer: object, smol: dict, force: bool = False, suppress: bool = False):
    """
    Remove a molecule from your molecule working set.

    Parameters
    ----------
    cmd_pointer: object
        The command pointer object.
    smol: dict
        The OpenAD molecule object to add.
    force: bool
        If True, add without confirming.
    suppress: bool
        If True, suppress success output.
    """

    if not smol:
        output_error("No molecule provided", return_val=False)
        return False

    # Name
    name = smol["name"]

    # Remove function.
    def _remove_mol():
        i = 0
        key, identifier_to_delete = get_best_available_identifier(smol)
        while cmd_pointer.molecule_list[i]["properties"][key] != identifier_to_delete:
            i = i + 1
        cmd_pointer.molecule_list.pop(i)
        if suppress is False:
            output_success(f"Molecule <yellow>{name}</yellow> was removed", pad=0, return_val=False)
        return True

    # Remove without confirming.
    if force:
        return _remove_mol()

    # Confirm before removing.
    smiles = get_best_available_smiles(smol)
    smiles_str = f" <reset>{smiles}</reset>" if smiles else ""
    if confirm_prompt(f"Remove molecule <green>{name}</green>{smiles_str} from your molecules working list?"):
        return _remove_mol()
    else:
        output_error(f"Molecule <yellow>{name}</yellow> was not removed", pad=0, return_val=False)
        return False


# endregion


# @@TODO: update merge_molecule command to use this function.
def merge_mols(openad_mol_1, openad_mol_2):
    """
    Merge two molecule objects into one.

    This adds parameters from 2nd molecule into the 1st, while
    preventing None values from overwriting existing values.
    """

    # Loop through all properties of the 2nd molecule.
    for key, value in openad_mol_2.items():
        # Skip if None
        if value is None:
            continue

        # If the property is a dictionary, merge it.
        if isinstance(value, dict):
            if key in openad_mol_1:
                openad_mol_1[key].update(value)
            else:
                openad_mol_1[key] = value

        # If the property is a list, merge it.
        elif isinstance(value, list):
            if key in openad_mol_1:
                # Temporary fix for synonyms, which are stored in a nested list in v1 dict.
                if key == "synonyms" and "Synonym" in openad_mol_1[key]:
                    openad_mol_1[key]["Synonym"].extend(value)
                else:
                    openad_mol_1[key].extend(value)
            else:
                openad_mol_1[key] = value

        # If the property is a string, merge it.
        elif isinstance(value, str):
            openad_mol_1[key] = value

    # Return
    return openad_mol_1


# @@TODO: merge function with merge_mols
def merge_molecule_properties(molecule_dict: dict, smol: dict):
    """
    Merge a dictionary with properties into a molecule's properties.

    Parameters
    ----------
    molecule_dict: dict
        The dictionary with properties to merge.
    smol: dict
        The OpenAD small molecule dictionary into which we merge the properties.
    """

    if smol is None:
        return None
    if "ROMol" in molecule_dict:
        del molecule_dict["ROMol"]
    if "subject" in molecule_dict:
        del molecule_dict["subject"]

    for key in molecule_dict:
        smol["properties"][key] = molecule_dict[key]
        smol["property_sources"][key] = {"source": "unknown", "date": pretty_date()}
        if key not in SMOL_PROPERTIES:
            SMOL_PROPERTIES.append(key)

    return smol


# @@TODO: merge function with merge_mols
def merge_molecule_REPLACE(merge_mol, mol):
    """merges a molecules property with those from a dictionary"""
    if mol is None:
        return None

    for key in merge_mol["properties"]:
        if key not in mol["properties"]:
            mol["properties"][key] = merge_mol["properties"][key]
            mol["property_sources"][key] = merge_mol["properties"][key]
        elif mol["properties"][key] is None:
            mol["properties"][key] = merge_mol["properties"][key]
            mol["property_sources"][key] = merge_mol["properties"][key]

    for x in merge_mol["analysis"]:
        if x not in mol["anaylsis"]:
            mol["anaylsis"].append()


############################################################
# region - DEPRECATED


# def molformat_v2(mol):
#     """
#     Create a slightly modified "v2" OpenAD molecule data format,
#     where identifiers are stored separately from properties. This
#     is how the GUI consumes molecules, and should be used elsewhere
#     in the application going forward as well. Please read comment
#     above new_molecule() for more information.
#     - - -
#     What it does:
#      1. Separate identifiers from properties.
#      2. Flatten the mol["synonyms"]["Synonym"] to mol["synonyms"]
#     """
#     if mol is None:
#         return None
#     mol_v2 = {}
#     mol_v2["identifiers"] = _get_identifiers(mol)
#     mol_v2["properties"] = deepcopy(mol.get("properties"))

#     # For the messy double-level synonyms key in v1 format.
#     if "Synonym" in mol.get("synonyms", {}):
#         synonyms = deepcopy(mol.get("synonyms", {}).get("Synonym", []))

#     # For cases like an SDF or CSV where everything is just stored flat.
#     elif "synonyms" in mol_v2["properties"]:
#         synonyms = mol_v2["properties"]["synonyms"]
#         del mol_v2["properties"]["synonyms"]

#     # For other cases, usually when no synonyms are present
#     # like when creating a fresh molecule from an identifier.
#     else:
#         synonyms = deepcopy(mol.get("synonyms", []))

#     # Synonyms may be stored as a string array (mdl/sdf), or a string with linebreaks (csv).
#     # print("\n\n* synonyms BEFORE", type(synonyms), "\n", synonyms)
#     if isinstance(synonyms, str):
#         try:
#             synonyms = ast.literal_eval(synonyms)
#         except (ValueError, SyntaxError):
#             synonyms = synonyms.split("\n") if "\n" in synonyms else [synonyms]
#     # print("\n\n* synonyms AFTER", type(synonyms), "\n", synonyms)

#     # Find name if it's missing.
#     if not mol_v2["identifiers"]["name"]:
#         mol_v2["identifiers"]["name"] = synonyms[0] if synonyms else None

#     mol_v2["synonyms"] = synonyms
#     mol_v2["analysis"] = deepcopy(mol.get("analysis"))
#     mol_v2["property_sources"] = deepcopy(mol.get("property_sources"))
#     mol_v2["enriched"] = deepcopy(mol.get("enriched"))

#     # # This removes all properties that are not in MOL_PROPERTIES
#     # # I'm not sure what the benefit is. This shouldn't be limited.
#     # mol_organized["properties"] = get_human_properties(mol)

#     # if "DS_URL" in mol_organized["properties"]:
#     #     mol_organized["properties"]["DS_URL"] = ""

#     # Remove identifiers from properties.
#     # - - -
#     # Create a lowercase version of the properties dictionary
#     # so we can scan for properties in a case-insensitive way.
#     molIdfrs = {k.lower(): v for k, v in mol_v2["identifiers"].items()}
#     for prop in list(mol_v2["properties"]):
#         if prop.lower() in molIdfrs:
#             del mol_v2["properties"][prop]

#     return mol_v2


# # Previous version of flatten_smol()
# def molformat_v2_to_v1(mol):
#     """
#     Convert the v2 OpenAD molecule object format back to the old v1 format.
#     """
#     if mol is None:
#         return None

#     mol_v1 = {}
#     mol_v1["name"] = deepcopy(mol.get("identifiers").get("name"))
#     mol_v1["properties"] = {**mol.get("identifiers"), **mol.get("properties")}
#     mol_v1["synonyms"] = {"Synonym": deepcopy(mol.get("synonyms"))}
#     mol_v1["analysis"] = deepcopy(mol.get("analysis"))
#     mol_v1["property_sources"] = deepcopy(mol.get("property_sources"))
#     mol_v1["enriched"] = deepcopy(mol.get("enriched"))
#     return mol_v1


# Unused. This adds an indice when it's missing, but there's no usecase
# other than dev-legacy example molsets that are missing an index.
def index_molset_file_async(path_absolute):
    """
    Add an index to each molecule of a molset file,
    without blocking the main thread.

    This is used to index a cached working copy of a molset
    right after it's created.

    Parameters
    ----------
    cache_path: str
        The path to the cached working copy of a molset.
    """

    async def _index_molset_file(cache_path):
        # Read
        async with aiofiles.open(cache_path, "r", encoding="utf-8") as f:
            content = await f.read()
        molset = json.loads(content)
        for i, mol in enumerate(molset):
            mol["index"] = i + 1
            molset[i] = mol
        # Write
        async with aiofiles.open(cache_path, "w", encoding="utf-8") as f:
            await f.write(json.dumps(molset, ensure_ascii=False, indent=4, cls=DecimalEncoder))

    asyncio.run(_index_molset_file(path_absolute))


# endregion
