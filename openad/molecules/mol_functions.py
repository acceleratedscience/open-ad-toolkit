"""
Molecule management functions
"""

import os
import time
import json
import copy
import shutil
import pandas
import asyncio
import aiofiles
import pubchempy as pcy
from datetime import datetime
from rdkit import Chem, rdBase
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.Descriptors import MolWt, ExactMolWt

from openad.helpers.output_msgs import msg
from openad.helpers.output import output_error, output_warning, output_success
from openad.helpers.files import open_file
from openad.helpers.general import confirm_prompt
from openad.helpers.spinner import spinner
from openad.helpers.json_decimal_encoder import DecimalEncoder

# The base for our molecule dictionary.
OPENAD_MOL_DICT = {
    "name": None,
    "synonyms": {},
    "properties": {},
    "property_sources": {},
    "sources": {},
    "commments": {},
    "analysis": [],
    "enriched": False,
}
MOL_NAME_INDEX = "name"
MOL_SMILES_INDEX = "smiles"
MOL_INCHI_INDEX = "inchi"
MOL_INCHIKEY_INDEX = "inchikey"
MOL_CID_INDEX = "cid"
MOL_SDF_INDEX = "sdf"
MOL_FORMULA = "formula"
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
    # "formal charge": "charge", # Duplicate but inconsistent
    "SOL_classification": "sol_classification",
    "SOL": "sol",
}
PROPERTY_SOURCES = {
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
    # "Mass-Exact": "molecular_weight", # ?
    "Weight-MonoIsotopic": "monoisotopic_mass",
    "Molecular Formula": "molecular_formula",
    "Topological-Polar Surface Area": "tpsa",
}
MOL_PROPERTIES = sorted(
    [
        "cid",
        "molecular_formula",
        "molecular_weight",
        "molecular_weight_exact",
        "canonical_smiles",
        "isomeric_smiles",
        "inchi",
        "inchikey",
        "iupac_name",
        "xlogp",
        "exact_mass",
        "monoisotopic_mass",
        "multipoles_3d",
        "tpsa",
        "complexity",
        "charge",
        "h_bond_donor_count",
        "h_bond_acceptor_count",
        "rotatable_bond_count",
        "heavy_atom_count",
        "isotope_atom_count",
        "atom_stereo_count",
        "defined_atom_stereo_count",
        "undefined_atom_stereo_count",
        "bond_stereo_count",
        "defined_bond_stereo_count",
        "undefined_bond_stereo_count",
        "covalent_unit_count",
        "volume_3d",
        "conformer_rmsd_3d",
        "conformer_model_rmsd_3d",
        "x_steric_quadrupole_3d",
        "y_steric_quadrupole_3d",
        "z_steric_quadrupole_3d",
        "feature_count_3d",
        "feature_acceptor_count_3d",
        "feature_donor_count_3d",
        "feature_anion_count_3d",
        "feature_cation_count_3d",
        "feature_ring_count_3d",
        "feature_hydrophobe_count_3d",
        "effective_rotor_count_3d",
        "conformer_count_3d",
        "pharmacophore_features_3d",
        "conformer_id_3d",
        "coordinate_type",
        "mmff94_energy_3d",
        "mmff94_partial_charges_3d",
        "multipoles_3d",
        "pharmacophore_features_3d",
        "sol_classification",
        "sol",
    ]
)

blocker = rdBase.BlockLogs()  # pylint: disable=c-extension-no-member
mol_name_cache = {}


def merge_molecule_properties(molecule_dict, mol):
    """merges a molecules property with those from a dictionary"""
    date_time = datetime.fromtimestamp(time.time())
    str_date_time = date_time.strftime("%d-%m-%Y, %H:%M:%S")
    if mol is None:
        return None
    if "ROMol" in molecule_dict:
        del molecule_dict["ROMol"]
    if "subject" in molecule_dict:
        del molecule_dict["subject"]

    for key in molecule_dict:
        mol["properties"][key] = molecule_dict[key]
        mol["property_sources"][key] = {"software": "custom", "date": str_date_time}
        if key not in MOL_PROPERTIES:
            MOL_PROPERTIES.append(key)
    return mol


def valid_smiles(input_molecule) -> bool:
    """Check if an string is valid SMILES definition."""

    # blocker = rdBase.BlockLogs()  # pylint: disable=c-extension-no-member
    try:
        m = Chem.MolFromSmiles(input_molecule, sanitize=False)  # pylint: disable=no-member
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


def valid_inchi(input_molecule) -> bool:
    """Check if a string is valid InChI molecule."""

    # blocker = rdBase.BlockLogs()  # pylint: disable=c-extension-no-member
    try:
        m = Chem.inchi.InchiToInchiKey(input_molecule)
    except Exception:  # pylint: disable=broad-exception-caught
        return False
    if m is None:
        return False
    else:
        return True


# TRASH - renamed to avoid lint errors with canonical_smiles as variable name in function
# def canonical_smiles(input_molecule):
#     """returns a cannonical smiles string based on rdkit"""
#     return Chem.MolToSmiles(Chem.MolFromSmiles(input_molecule), True)  # pylint: disable=no-member


def canonicalize(input_smiles):
    """
    Turn any SMILES into its canonical equivalent per RDKit.
    """
    return Chem.MolToSmiles(Chem.MolFromSmiles(input_smiles), isomericSmiles=True)  # pylint: disable=no-member


def retrieve_mol_from_list(cmd_pointer, identifier, ignore_synonyms=False):
    """retrieves a molecule from the working list"""

    openad_mol = find_mol_in_list(identifier, cmd_pointer.molecule_list, ignore_synonyms=ignore_synonyms)
    if openad_mol is not None:
        return openad_mol.copy()
    return None

    # TRASH
    # for mol in cmd_pointer.molecule_list:
    #     m = is_molecule(mol, molecule)

    #     if m is not None:
    #         return m.copy()


# TODO: very convoluted way to call _get_mol()
def retrieve_mol(molecule):
    """Fetch molecule from PubChem"""
    success, mol, comp = get_mol_from_name(molecule)
    if success:
        return mol

    success, mol, comp = get_mol_from_inchi(molecule)
    if success:
        return mol

    success, mol, comp = get_mol_from_smiles(molecule)
    if success:
        return mol

    # Commented out until getting no time outs from pubchem.
    # success, mol, comp = get_mol_from_formula(molecule)
    # if success:
    #    return mol

    success, mol, comp = get_mol_from_inchikey(molecule)
    if success:
        return mol

    success, mol, comp = get_mol_from_cid(molecule)
    if success:
        return mol
    return None


def get_mol_from_inchi(inchi_str: str):
    """return pubchem molecule data based on inchi"""

    if Chem.MolFromInchi(inchi_str):
        return _get_mol(inchi_str, MOL_INCHI_INDEX)
    else:
        return False, "Invalid Inchi", None


def get_mol_from_inchikey(inchikey_str: str):
    """return pubchem molecule data based on inchi key"""
    return _get_mol(inchikey_str, MOL_INCHIKEY_INDEX)


def get_mol_from_smiles(smiles_str: str):
    """return pubchem molecule data based on smiles"""

    if valid_smiles(smiles_str):
        # print("getting smiles")

        success, openad_mol, molecule = _get_mol(smiles_str, MOL_SMILES_INDEX)
        if openad_mol is not None:
            openad_mol["properties"]["canonical_smiles"] = Chem.MolToSmiles(Chem.MolFromSmiles(smiles_str))

        return success, openad_mol, molecule
    else:
        return False, "Invalid Smiles", None


def get_mol_from_name(mol_name: str):
    """return pubchem molecule data based on molecule name"""
    return _get_mol(mol_name, MOL_NAME_INDEX)


def get_mol_from_formula(mol_formula: str):
    """return pubchem molecule data based on molecule formula"""
    return _get_mol(mol_formula, MOL_FORMULA)


def get_mol_from_cid(mol_cid: str):
    """return pubchem molecule data based on molecule cid"""
    return _get_mol(mol_cid, MOL_CID_INDEX)


def _get_mol(mol_id, mol_id_type):
    """gets molecule based on provided identifier"""
    cid = None
    try:
        if mol_id_type == MOL_CID_INDEX:
            if not mol_id.isdigit():
                return False, None, None
            cid = int(mol_id)

        else:
            compounds = pcy.get_compounds(mol_id, mol_id_type)
            if len(compounds) == 0:
                return False, None, None
            else:
                cid = compounds[0].cid
            if cid is None:
                return False, None, None
        mol_pcy = pcy.Compound.from_cid(cid).to_dict()
        openad_mol = copy.deepcopy(OPENAD_MOL_DICT)

        if mol_pcy:
            if mol_id_type == MOL_NAME_INDEX:
                openad_mol["name"] = mol_id
            else:
                openad_mol["name"] = mol_pcy["iupac_name"]
            openad_mol["properties"] = mol_pcy
            openad_mol["sources"]["pubchem"] = mol_pcy

            for x in MOL_PROPERTIES:
                openad_mol["property_sources"][x] = {"source": "pubchem"}

                for key, value in PROPERTY_SOURCES.items():
                    if value == x:
                        if len(key.split("-")) > 0:
                            for y in openad_mol["sources"]["pubchem"]["record"]["props"]:
                                if "label" not in y["urn"]:
                                    pass
                                elif y["urn"]["label"] == key.split("-")[0] and "name" not in y["urn"]:
                                    openad_mol["property_sources"][x] = y["urn"]
                                elif y["urn"]["label"] == key.split("-")[0] and y["urn"]["name"] == key.split("-")[1]:
                                    openad_mol["property_sources"][x] = y["urn"]
            names = pcy.get_synonyms(openad_mol["name"], "name")
            if len(names) > 0:
                openad_mol["synonyms"] = names[0]
            openad_mol["enriched"] = True

            return True, openad_mol, mol_pcy
    except Exception as e:
        print(e)
        return False, None, None


def get_properties(mol):
    """pulls properties from list"""
    properties_dict = {}
    for mol_property in MOL_PROPERTIES:
        if mol_property in mol["properties"]:
            properties_dict[mol_property] = mol["properties"][mol_property]
    return properties_dict


def get_identifiers(mol):
    """pulls the identifiers from a molecule"""
    identifier_dict = {}

    # Create a lowercase version of the properties dictionary
    # so we can scan for properties in a case-insensitive way.
    molProps = {k.lower(): v for k, v in mol["properties"].items()}

    # SDF files will have 'name' as a property, while our old mol format has it as a main key.
    identifier_dict["name"] = mol.get("name") or molProps.get("name")
    identifier_dict["inchi"] = molProps.get("inchi")
    identifier_dict["inchikey"] = molProps.get("inchikey")
    identifier_dict["canonical_smiles"] = molProps.get("canonical_smiles")
    identifier_dict["isomeric_smiles"] = molProps.get("isomeric_smiles")
    identifier_dict["molecular_formula"] = molProps.get("molecular_formula")
    identifier_dict["cid"] = molProps.get("cid")

    # We don't use the unspecified "smiles" property,
    # but when parsing an SDF file, this may be a key.
    identifier_dict["smiles"] = molProps.get("smiles")

    return identifier_dict


def read_molset_from_cache(cmd_pointer, cache_id):
    """
    Read a cached molset file from disk.

    Parameters
    ----------
    cmd_pointer: object
        The command pointer object, used to fetch the workspace path.
    cache_id: str
        The cache ID of the molset file.
    """

    # Read file from cache.
    cache_path = assemble_cache_path(cmd_pointer, "molset", cache_id)
    molset, err_code = get_molset_mols(cache_path)

    # Return error
    if err_code:
        raise ValueError(f"Failed to read molset {cache_id} from cache")
    else:
        return molset


def get_molset_mols(path_absolute):
    """
    Return the list of molecules from a molset file,
    with an index added to each molecule.

    Parameters
    ----------
    path_absolute: str
        The absolute path to the molset file.

    Returns: data, err_code
    """
    # Read file contents
    molset, err_code = open_file(path_absolute, return_err="code")

    return molset, err_code


def molformat_v2(mol):
    """
    Create a slightly modified "v2" OpenAD molecule data format,
    where identifiers are stored separately from properties. This
    is how the GUI consumes molecules, and should be used elsewhere
    in the application going forward as well. Please read comment
    above new_molecule() for more information.
    - - -
    What it does:
     1. Separate identifiers from properties.
     2. Flatten the mol["synonyms"]["Synonym"] to mol["synonyms"]
    """
    if mol is None:
        return None
    mol_v2 = {}
    mol_v2["identifiers"] = get_identifiers(mol)
    mol_v2["properties"] = copy.deepcopy(mol.get("properties"))
    # For the messy double-level synonyms key in v1 format.
    mol_v2["synonyms"] = copy.deepcopy(mol.get("synonyms", {}).get("Synonym", []))
    # For other cases like an SDF or CSV where everything is just stored flat.
    if not mol_v2["synonyms"] and "synonyms" in mol_v2["properties"]:
        mol_v2["synonyms"] = str(mol_v2["properties"]["synonyms"]).split("\n")
        del mol_v2["properties"]["synonyms"]
    mol_v2["analysis"] = copy.deepcopy(mol.get("analysis"))
    mol_v2["property_sources"] = copy.deepcopy(mol.get("property_sources"))
    mol_v2["enriched"] = copy.deepcopy(mol.get("enriched"))

    # # This removes all properties that are not in MOL_PROPERTIES
    # # I'm not sure what the benefit is. This shouldn't be limited.
    # mol_organized["properties"] = get_properties(mol)

    # if "DS_URL" in mol_organized["properties"]:
    #     mol_organized["properties"]["DS_URL"] = ""

    # Remove identifiers from properties.
    # - - -
    # Create a lowercase version of the properties dictionary
    # so we can scan for properties in a case-insensitive way.
    molIdfrs = {k.lower(): v for k, v in mol_v2["identifiers"].items()}
    for prop in list(mol_v2["properties"]):
        if prop.lower() in molIdfrs:
            del mol_v2["properties"][prop]

    return mol_v2


def molformat_v2_to_v1(mol):
    """
    Convert the v2 OpenAD molecule object format back to the old v1 format.
    """
    if mol is None:
        return None
    mol_v1 = copy.deepcopy(OPENAD_MOL_DICT)
    mol_v1["name"] = copy.deepcopy(mol.get("identifiers").get("name"))
    mol_v1["properties"] = {**mol.get("identifiers"), **mol.get("properties")}
    mol_v1["synonyms"]["Synonym"] = copy.deepcopy(mol.get("synonyms"))
    mol_v1["analysis"] = copy.deepcopy(mol.get("analysis"))
    mol_v1["property_sources"] = copy.deepcopy(mol.get("property_sources"))
    mol_v1["enriched"] = copy.deepcopy(mol.get("enriched"))
    return mol_v1


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
def new_molecule(inchi_or_smiles: str = None, mol_rdkit: Mol = None, name: str = None):
    """
    Create a basic molecule object without relying on API calls.

    Parameters
    ----------
    inchi_or_smiles: str
        Source option A: An InChI or SMILES identifier.
    mol_rdkit: RDKit molecule object
        Source option B: An RDKit molecule object.
    name: str
        Optional name for the molecule.
    """

    mol = copy.deepcopy(OPENAD_MOL_DICT)
    date_time = datetime.fromtimestamp(time.time())
    str_date_time = date_time.strftime("%d-%m-%Y, %H:%M:%S")

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

    # This is creating a bunch of unnecessary empty fields... marked for deletion
    # # Create empty property fields
    # for prop in MOL_PROPERTIES:
    #     mol["properties"][prop] = None

    # Store identifiers
    mol["name"] = name

    mol["properties"]["inchi"] = Chem.MolToInchi(mol_rdkit)  # Alt: Chem.rdinchi.MolToInchi(mol_rdkit)[0]
    mol["property_sources"]["inchi"] = {"software": "rdkit", "date": str_date_time}

    mol["properties"]["inchikey"] = Chem.inchi.InchiToInchiKey(mol["properties"]["inchi"])
    mol["property_sources"]["inchikey"] = {"software": "rdkit", "date": str_date_time}

    mol["properties"]["canonical_smiles"] = Chem.MolToSmiles(mol_rdkit)  # pylint: disable=no-member
    mol["property_sources"]["canonical_smiles"] = {"software": "rdkit", "date": str_date_time}

    mol["properties"]["isomeric_smiles"] = Chem.MolToSmiles(mol_rdkit, isomericSmiles=True)  # pylint: disable=no-member
    mol["property_sources"]["isomeric_smiles"] = {"software": "rdkit", "date": str_date_time}

    mol["properties"]["molecular_formula"] = Chem.rdMolDescriptors.CalcMolFormula(mol_rdkit)
    mol["property_sources"]["molecular_formula"] = {"software": "rdkit", "date": str_date_time}

    mol["properties"]["molecular_weight"] = MolWt(mol_rdkit)
    mol["property_sources"]["molecular_weight"] = {"software": "rdkit", "date": str_date_time}

    # Disabled this for consistency, because molecular_weight_exact is not included in PubChem data.
    # See: pcy.Compound.from_cid(cid).to_dict()
    # mol["properties"]["molecular_weight_exact"] = ExactMolWt(mol_rdkit)
    # mol["property_sources"]["molecular_weight_exact"] = {"software": "rdkit", "date": str_date_time}

    # So the UI can recognize when a molecule has been enriched.
    mol["enriched"] = False

    return mol


def create_molset_cache_file(cmd_pointer, molset=None, path_absolute=None):
    """
    Store molset as a cached file, so we can manipulate it in the GUI.

    Returns cache_id
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


def assemble_cache_path(cmd_pointer, file_type, cache_id):
    """
    Compile the file path to a cached working copy of a file.

    Parameters:
    -----------
    cmd_pointer: object
        The command pointer object, used to fetch the workspace path.
    file_type: 'molset'
        The type of file, used to name the cache file. For now only molset.
    """
    workspace_path = cmd_pointer.workspace_path()
    return f"{workspace_path}/._openad/wc_cache/{file_type}-{cache_id}.json"


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


def df_has_molecules(df):
    """
    Check if a dataframe has molecules.
    """
    colsLowercase = [col.lower() for col in df.columns]
    if "smiles" in colsLowercase:
        return True
    elif "inchi" in colsLowercase:
        return True
    else:
        return False


def mol_from_identifier(cmd_pointer, identifier, mol_name=None, basic=False):
    """
    Fetch a OpenAD molecule dict from an identifier.
    Either from memory or PubChem.

    Parameters
    ----------
    cmd_pointer: object
        The command pointer object.
    identifier: str
        The molecule identifier to search for.
        Valid inputs: name, CID, InChI, InChIKey, SMILES.
    basic: bool
        If True, create a basic molecule dict with RDKit (no API calls).
    mol_name: str
        Optional name for the molecule.
    """
    # Create basic molecule dict with RDKit (no API calls).
    if basic is True:
        openad_mol = new_molecule(identifier, name=mol_name)

    # Create enriched molecule dict.
    else:
        # Fetch from memory.
        if (
            cmd_pointer.last_external_molecule is not None
            and find_mol_in_list(identifier, [cmd_pointer.last_external_molecule]) is not None
        ):
            openad_mol = cmd_pointer.last_external_molecule

        # Fetch from PubChem.
        else:
            openad_mol = retrieve_mol(identifier)  # TODO: very convoluted way to call _get_mol()
            if openad_mol is None:
                openad_mol = new_molecule(identifier, name=identifier)

    # Fail - invalid.
    if openad_mol is None:
        output_error("Unable to identify molecule", return_val=False)
        if basic is True:
            output_error("Valid SMILES string required", return_val=False)

    return openad_mol


def mymols_add(cmd_pointer, openad_mol, force=False, suppress=False):
    """
    Add a molecule to your molecule bookmarks.

    Parameters
    ----------
    cmd_pointer: object
        The command pointer object.
    openad_mol: dict
        The OpenAD molecule object to add.
    force: bool
        If True, add without confirming.
    suppress: bool
        If True, suppress success output.
    """

    if not openad_mol:
        output_error("No molecule provided", return_val=False)
        return False

    # Fail - already in list.
    if retrieve_mol_from_list(cmd_pointer, openad_mol["properties"]["canonical_smiles"]) is not None:
        output_error("Molecule already in list: " + openad_mol["properties"]["canonical_smiles"], return_val=False)
        return True

    # Name
    name = openad_mol["name"]

    # Add function
    def _add_mol():
        cmd_pointer.molecule_list.append(openad_mol.copy())
        if suppress is False:
            output_success(f"Molecule <yellow>{name}</yellow> was added", pad=0, return_val=False)
        return True

    # Add without confirming.
    if force:
        return _add_mol()

    # Confirm before adding.
    smiles = get_best_available_smiles(openad_mol)
    smiles_str = f" <reset>{smiles}</reset>" if smiles else ""
    if confirm_prompt(f"Add molecule <green>{name}</green>{smiles_str} to your molecules working list?"):
        return _add_mol()
    else:
        output_error(f"Molecule <yellow>{name}</yellow> was not added", pad=0, return_val=False)
        return False


def mymols_remove(cmd_pointer, openad_mol, force=False, suppress=False):
    """
    Remove a molecule from your molecule bookmarks.

    Parameters
    ----------
    cmd_pointer: object
        The command pointer object.
    openad_mol: dict
        The OpenAD molecule object to add.
    force: bool
        If True, add without confirming.
    suppress: bool
        If True, suppress success output.
    """

    if not openad_mol:
        output_error("No molecule provided", return_val=False)
        return False

    # Name
    name = openad_mol["name"]

    # Remove function.
    def _remove_mol():
        i = 0
        key, identifier_to_delete = get_best_available_identifier(openad_mol)
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
    smiles = get_best_available_smiles(openad_mol)
    smiles_str = f" <reset>{smiles}</reset>" if smiles else ""
    if confirm_prompt(f"Remove molecule <green>{name}</green>{smiles_str} from your molecules working list?"):
        return _remove_mol()
    else:
        output_error(f"Molecule <yellow>{name}</yellow> was not removed", pad=0, return_val=False)
        return False


def find_mol_in_list(identifier, molset, ignore_synonyms=False):
    """
    Scan a molset for a given identifier.
    Returns the molecule dict if found, otherwise None.

    Parameters
    ----------
    identifier: str
        The molecule identifier to search for.
        Valid inputs: name, CID, InChI, InChIKey, SMILES.
    molset: list
        A list of OpenAD molecule objects.
    ignore_synonyms: bool
        If True, ignore synonyms in the search.
        This is only used when renaming a molecule to one of its synonyms.
        Without it, the search would return the original molecule and abort.

    """

    for openad_mol in molset:
        # To support both v1 and v2 formats (see molformat_v2).
        identifiers_dict = openad_mol.get("identifiers")  # v2
        if not identifiers_dict:
            identifiers_dict = openad_mol.get("properties")  # v1
        if not identifiers_dict:
            output_error("find_mol_in_list(): Invalid molset input", return_val=False)
            return None
        synonyms = openad_mol.get("synonyms", {}).get("Synonym") or openad_mol.get("synonyms", [])

        # Name match v1
        if identifier.upper() == (openad_mol.get("name", "") or "").upper():
            return openad_mol

        # Name match v2 format
        if identifier.upper() == identifiers_dict.get("name", "").upper():
            return openad_mol

        # CID match
        try:
            if int(identifier) == int(identifiers_dict.get("cid")):
                return openad_mol
        except Exception:  # pylint: disable=broad-except
            pass

        # InChI match
        if identifier == identifiers_dict.get("inchi"):
            return openad_mol

        # InChIKey match
        if identifier == identifiers_dict.get("inchikey"):
            return openad_mol

        # SMILES match - isomeric
        if (
            identifiers_dict.get("isomeric_smiles") is not None
            and identifier.upper() == identifiers_dict.get("isomeric_smiles", "").upper()
        ):
            return openad_mol

        # SMILES match - canonical
        try:
            if canonicalize(identifier) == canonicalize(identifiers_dict.get("canonical_smiles")):
                return openad_mol
        except Exception:  # pylint: disable=broad-except
            pass

        # Synonym match
        if not ignore_synonyms:
            for syn in synonyms:
                # print(syn, identifier)
                if identifier.upper() == syn.upper():
                    return openad_mol

    # Fail
    return None


def get_best_available_identifier(openad_mol):  # v1 or v2 molecule object
    """
    Get whatever identifier is available from a molecule.
    """

    identifiers_dict = openad_mol.get("identifiers")  # For v2 molecule objects.
    if not identifiers_dict:
        identifiers_dict = openad_mol.get("properties")  # For v1 molecule objects.

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


def get_best_available_smiles(openad_mol):  # v1 or v2 molecule object
    """
    Get the best available SMILES string from a molecule.
    """

    identifiers_dict = openad_mol.get("identifiers")  # For v2 molecule objects.
    if not identifiers_dict:
        identifiers_dict = openad_mol.get("properties")  # For v1 molecule objects.

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


# Todo: update merge_molecule command to use this function.
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
                openad_mol_1[key].extend(value)
            else:
                openad_mol_1[key] = value

        # If the property is a string, merge it.
        elif isinstance(value, str):
            openad_mol_1[key] = value

    # Return
    return openad_mol_1


def normalize_mol_df(mol_df: pandas.DataFrame, batch=False):
    """
    Normalize the column names of a molecule dataframe
    """
    has_name = False
    contains_name = None

    for i in mol_df.columns:
        # Find the name column.
        if str(i.upper()) == "NAME" or str(i.lower()) == "chemical_name":
            has_name = True
        if contains_name is None and "NAME" in str(i.upper()):
            contains_name = i
        if contains_name is None and "CHEMICAL_NAME" in str(i.upper()):
            contains_name = i

        # Normalize any columns we'll be referring to later.
        if str(i.upper()) == "SMILES":
            mol_df.rename(columns={i: "SMILES"}, inplace=True)
        if str(i.upper()) == "ROMOL":
            mol_df.rename(columns={i: "ROMol"}, inplace=True)
        if str(i.upper()) == "IMG":
            mol_df.rename(columns={i: "IMG"}, inplace=True)
        if i in INPUT_MAPPINGS:
            mol_df.rename(columns={i: INPUT_MAPPINGS[i]}, inplace=True)

    # Normalize name column.
    if has_name == False and contains_name is not None:
        mol_df.rename(columns={contains_name: "NAME"}, inplace=True)

    # Add names when missing.
    try:
        if has_name is False:
            output_warning(msg("no_m2g_name_column"))

            if not batch:
                spinner.start("Downloading names")

            mol_df["NAME"] = "unknown"
            for i in mol_df.itertuples():
                mol_df.loc[i.Index, "NAME"] = _smiles_to_iupac(mol_df.loc[i.Index, "SMILES"])

            if not batch:
                spinner.succeed("Names downloaded")
                spinner.start()
                spinner.stop()

    except BaseException as err:
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

    import pubchempy

    if smiles in mol_name_cache:
        return mol_name_cache[smiles]
    try:
        compounds = pubchempy.get_compounds(smiles, namespace="smiles")
        match = compounds[0]
        mol_name_cache[smiles] = str(match)
    except BaseException:
        match = smiles
    return str(match)
