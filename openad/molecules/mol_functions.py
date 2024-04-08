"""Implements molecule management functions"""

import time
import copy
import pubchempy as pcy
from datetime import datetime
from rdkit import Chem, rdBase
from rdkit.Chem.Descriptors import MolWt, ExactMolWt
from openad.helpers.output import output_text, output_table, output_warning, output_error
from openad.helpers.files import open_file

blocker = rdBase.BlockLogs()

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
    "molecular weight": "molecular_weight",
    "monoisotopic mass": "monoisotopic_mass",
    "topological polar surface area": "tpsa",
    "heavy atom count": "heavy_atom_count",
    "formal charge": "formal_change",
    "isotope atom count": "isotope_atom_count",
    "defined atom stereocenter count": "defined_atom_stereo_count",
    "undefined atom stereocenter count": "undefined_atom_stereo_count",
    "covalently-bonded unit count": "covalent_unit_count",
    "compound is canonicalized": "compound_canonicalized",
    "formal charge": "charge",
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
    "Mass-Exact": "molecular_weight",
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


def merge_molecule_properties(molecule_dict, mol):
    """merges a molecules property with those from a dictionary"""
    date_time = datetime.fromtimestamp(time.time())
    str_date_time = date_time.strftime("%d-%m-%Y, %H:%M:%S")
    if mol is None:
        return None
    if "ROMol" in molecule_dict:
        del molecule_dict["ROMol"]

    for key in molecule_dict:
        mol["properties"][key] = molecule_dict[key]
        mol["property_sources"][key] = {"software": "custom", "date": str_date_time}
        if key not in MOL_PROPERTIES:
            MOL_PROPERTIES.append(key)
    return mol


def valid_smiles(input_molecule) -> bool:
    """Check if an string is valid SMILES definition."""

    blocker = rdBase.BlockLogs()  # pylint: disable=c-extension-no-member
    try:
        m = Chem.MolFromSmiles(input_molecule, sanitize=False)  # pylint: disable=no-member
    except:
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

    blocker = rdBase.BlockLogs()  # pylint: disable=c-extension-no-member
    try:
        m = Chem.inchi.InchiToInchiKey(input_molecule)
    except:
        return False
    if m is None:
        return False
    else:
        return True


def canonical_smiles(input_molecule):
    """returns a cannonical smiles string based on rdkit"""
    return Chem.MolToSmiles(Chem.MolFromSmiles(input_molecule), True)  # pylint: disable=no-member (false positive)


def get_mol_from_inchi(inchi_str: str):
    """return pubchem molecule data based on inchi"""

    if Chem.MolFromInchi(inchi_str):
        return get_mol(inchi_str, MOL_INCHI_INDEX)
    else:
        return False, "Invalid Inchi", None


def get_mol_from_inchikey(inchikey_str: str):
    """return pubchem molecule data based on inchi key"""
    return get_mol(inchikey_str, MOL_INCHIKEY_INDEX)


def get_mol_from_smiles(smiles_str: str):
    """return pubchem molecule data based on smiles"""
    if valid_smiles(smiles_str):
        # print("getting smiles")
        return get_mol(smiles_str, MOL_SMILES_INDEX)
    else:
        return False, "Invalid Smiles", None


def get_mol_from_name(mol_name: str):
    """return pubchem molecule data based on molecule name"""
    return get_mol(mol_name, MOL_NAME_INDEX)


def get_mol_from_formula(mol_formula: str):
    """return pubchem molecule data based on molecule formula"""
    return get_mol(mol_formula, MOL_FORMULA)


def get_mol_from_cid(mol_cid: str):
    """return pubchem molecule data based on molecule cid"""
    return get_mol(mol_cid, MOL_CID_INDEX)


def get_mol(mol_id, mol_id_type):
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
        molecule = pcy.Compound.from_cid(cid).to_dict()
        openad_mol = {
            "name": None,
            "synonyms": {},
            "properties": {},
            "property_sources": {},
            "sources": {},
            "commments": {},
            "analysis": [],
        }

        if molecule:
            if mol_id_type == MOL_NAME_INDEX:
                openad_mol["name"] = mol_id
            else:
                openad_mol["name"] = molecule["iupac_name"]
            openad_mol["properties"] = molecule
            openad_mol["sources"]["pubchem"] = molecule

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

            return True, openad_mol, molecule
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

    identifier_dict["name"] = mol.get("name")
    identifier_dict["inchi"] = molProps.get("inchi")
    identifier_dict["inchikey"] = molProps.get("inchikey")
    identifier_dict["canonical_smiles"] = molProps.get("canonical_smiles")
    identifier_dict["isomeric_smiles"] = molProps.get("isomeric_smiles")
    identifier_dict["formula"] = molProps.get("molecular_formula")
    identifier_dict["cid"] = molProps.get("cid")

    # We don't use the unspecified "smiles" property,
    # but when parsing an SDF file, this may be a key.
    identifier_dict["smiles"] = molProps.get("smiles")

    return identifier_dict


# This creates a slightly modified "v2" molecule data format, where
# identifiers are stored separately from properties. This is
# how the GUI consumes molecules, and should be used elsewhere
# in the application going forward as well. Please read comment
# above new_molecule() for more information.
# - - -
# What it does:
#  1. Separate identifiers from properties.
#  2. Flatten the mol["synonyms"]["Synonym"] to mol["synonyms"]
def molformat_v2(mol):
    if mol is None:
        return None
    mol_organized = {}
    mol_organized["identifiers"] = get_identifiers(mol)
    mol_organized["synonyms"] = copy.deepcopy(mol.get("synonyms", {}).get("Synonym", []))
    mol_organized["properties"] = copy.deepcopy(mol.get("properties"))
    mol_organized["analysis"] = copy.deepcopy(mol.get("analysis"))
    mol_organized["property_sources"] = copy.deepcopy(mol.get("property_sources"))
    mol_organized["enriched"] = copy.deepcopy(mol.get("enriched"))

    # # This removes all properties that are not in MOL_PROPERTIES
    # # I'm not sure what the benefit is. This shouldn't be limited.
    # mol_organized["properties"] = get_properties(mol)

    # if "DS_URL" in mol_organized["properties"]:
    #     mol_organized["properties"]["DS_URL"] = ""

    # Remove identifiers from properties.
    # - - -
    # Create a lowercase version of the properties dictionary
    # so we can scan for properties in a case-insensitive way.
    molIdfrs = {k.lower(): v for k, v in mol_organized["identifiers"].items()}
    for prop in list(mol_organized["properties"]):
        if prop.lower() in molIdfrs:
            del mol_organized["properties"][prop]

    return mol_organized


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
def new_molecule(inchi_or_smiles: str, name: str = None):
    """
    Create a basic molecule object without relying on API calls
    """

    mol = copy.deepcopy(OPENAD_MOL_DICT)
    date_time = datetime.fromtimestamp(time.time())
    str_date_time = date_time.strftime("%d-%m-%Y, %H:%M:%S")

    # Create RDKit molecule object
    try:
        mol_rdkit = Chem.MolFromInchi(inchi_or_smiles)
        if not mol_rdkit:
            mol_rdkit = Chem.MolFromSmiles(inchi_or_smiles)  # pylint: disable=no-member (false positive)
        if not mol_rdkit:
            mol_rdkit = Chem.MolFromInchi("InChI=1S/" + inchi_or_smiles)
        if not mol_rdkit:
            return None
    except Exception:  # pylint: disable=broad-exception-caught
        return None

    # Create empty property fields
    for prop in MOL_PROPERTIES:
        mol["properties"][prop] = None

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

    mol["properties"]["molecular_weight_exact"] = ExactMolWt(mol_rdkit)
    mol["property_sources"]["molecular_weight_exact"] = {"software": "rdkit", "date": str_date_time}

    # So the UI can recognize when a molecule has been enriched.
    mol["enriched"] = True

    return mol


# Create svg code from .
def mol2svg(mol_rdkit, highlight=None):
    if highlight:
        substructure = Chem.MolFromSmarts(highlight)  # pylint: disable=no-member
        matches = mol_rdkit.GetSubstructMatches(substructure)

        # Flatten the tuple of tuples into a list of atom indices
        highlight_atoms = [atom_index for match in matches for atom_index in match]
    else:
        highlight_atoms = None

    mol_drawer = Chem.Draw.MolDraw2DSVG(400, 300)  # pylint: disable=no-member
    mol_drawer.DrawMolecule(mol_rdkit, highlightAtoms=highlight_atoms)
    mol_drawer.FinishDrawing()
    return mol_drawer.GetDrawingText()


# Create sdf code.
def mol2sdf(mol_rdkit):
    # Add hydrogen atoms, which are displayed as spikes in the 3D viz.
    mol_rdkit = Chem.AddHs(mol_rdkit)  # pylint: disable=no-member
    # Generate 3D coordinates.
    Chem.rdDistGeom.EmbedMolecule(mol_rdkit)

    # Generate SDF format.
    mol_sdf = Chem.MolToMolBlock(mol_rdkit)  # pylint: disable=no-member
    return mol_sdf


# Not used, for testing
def mol2xyz(mol_rdkit):
    mol_xyz = Chem.rdmolfiles.MolToXYZBlock(mol_rdkit)
    return mol_xyz


# Not used, for testing
def mol2pdb(mol_rdkit):
    mol_pdb = Chem.rdmolfiles.MolToPDBBlock(mol_rdkit, flavor=32)
    return mol_pdb


# Check if a dataframe has molecules
def df_has_molecules(df):
    colsLowercase = [col.lower() for col in df.columns]
    if "smiles" in colsLowercase:
        return True
    elif "inchi" in colsLowercase:
        return True
    else:
        return False


if __name__ == "__main__":
    # get_mol_basic("InChI=1S/C13H18O2/c1-9(2)8-11-4-6-12(7-5-11)10(3)13(14)15/h4-7,9-10H,8H2,1-3H3,(H,14,15)")
    # get_mol_basic("ibuprofen")
    success, molec, cmp = get_mol_from_name("ibuprofen")
    print(success)
    if success:
        print(get_identifiers(molec))
        print(get_properties(molec))
