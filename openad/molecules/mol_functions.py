"""Implements molecule management functions"""

from datetime import datetime
import time
import pubchempy as pcy
from rdkit import Chem, rdBase
from openad.helpers.output import output_text, output_table, output_warning, output_error

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
}

MOL_NAME_INDEX = "name"
MOL_SMILES_INDEX = "smiles"
MOL_INCHI_INDEX = "inchi"
MOL_INCHIKEY_INDEX = "inchikey"
MOL_CID_INDEX = "cid"
MOL_SDF_INDEX = "sdf"
MOL_FORMULA = "formula"
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

MOL_PROPERTIES = [
    "cid",
    "molecular_formula",
    "molecular_weight",
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
]


def new_molecule(Name: str, smiles: str):
    """creates a new molecule object"""
    try:
        new_smiles = None
        rdkit_mol = None
        mol_weight = None
        inchi = None
        inchikey = None
        formula = None

        new_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(smiles), False)
        rdkit_mol = Chem.MolFromSmiles(new_smiles)
        formula = Chem.rdMolDescriptors.CalcMolFormula(rdkit_mol)
        inchi = Chem.rdinchi.MolToInchi(rdkit_mol)[0]
        inchikey = Chem.inchi.InchiToInchiKey(inchi)

        # mol_weight = Chem.Descriptors.ExactMolWt(rdkit_mol)

    except:
        if new_smiles is None:
            return None
    mol = {
        "name": Name,
        "synonyms": {},
        "property_sources": {},
        "properties": {},
        "commments": {},
        "analysis": [],
        "sources": [],
    }
    for i in MOL_PROPERTIES:
        mol["properties"][i] = None

    date_time = datetime.fromtimestamp(time.time())
    str_date_time = date_time.strftime("%d-%m-%Y, %H:%M:%S")
    mol["properties"]["molecular_weight"] = mol_weight
    mol["property_sources"]["molecular_weight"] = {"software": "rdkit", "date": str_date_time}
    mol["properties"]["inchi"] = inchi
    mol["property_sources"]["inchi"] = {"software": "rdkit", "date": str_date_time}
    mol["properties"]["inchikey"] = inchikey
    mol["property_sources"]["inchikey"] = {"software": "rdkit", "date": str_date_time}
    mol["properties"]["canonical_smiles"] = new_smiles
    mol["property_sources"]["canonical_smiles"] = {"software": "rdkit", "date": str_date_time}
    mol["properties"]["molecular_formula"] = formula
    mol["property_sources"]["molecular_formula"] = {"software": "rdkit", "date": str_date_time}
    return mol


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
    return Chem.MolToSmiles(Chem.MolFromSmiles(input_molecule), True)


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
        openad_mol = OPENAD_MOL_DICT.copy()

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

    identifier_dict["name"] = mol["name"]
    identifier_dict["cid"] = mol["properties"]["cid"]
    identifier_dict["inchi"] = mol["properties"]["inchi"]
    identifier_dict["inchikey"] = mol["properties"]["inchikey"]
    identifier_dict["isomeric_smiles"] = mol["properties"]["isomeric_smiles"]
    identifier_dict["canonical_smiles"] = mol["properties"]["canonical_smiles"]
    identifier_dict["formula"] = mol["properties"]["molecular_formula"]
    return identifier_dict


# Takes any identifier and creates a minimal molecule object,
# ideally without connecting to PubChem.
def get_mol_basic(molecule_identifier):
    mol = OPENAD_MOL_DICT.copy()

    # Fill in the basic identifiers
    mol_rdkit = Chem.MolFromInchi(molecule_identifier)
    if mol_rdkit:
        # Prepend 'InChI=' when missing.
        if not molecule_identifier.startswith("InChI="):
            molecule_identifier = "InChI=" + molecule_identifier

        mol["properties"]["inchi"] = molecule_identifier
        mol["properties"]["inchikey"] = Chem.inchi.InchiToInchiKey(molecule_identifier)
        mol["properties"]["canonical_smiles"] = Chem.MolToSmiles(mol_rdkit)
        mol["properties"]["isomeric_smiles"] = Chem.MolToSmiles(mol_rdkit, isomericSmiles=True)
        return mol
    else:
        mol_rdkit = Chem.MolFromSmiles(molecule_identifier)
    if mol_rdkit:
        mol["properties"]["inchi"] = Chem.MolToInchi(mol_rdkit)
        mol["properties"]["inchikey"] = Chem.inchi.InchiToInchiKey(mol["properties"]["inchi"])
        mol["properties"]["canonical_smiles"] = Chem.MolToSmiles(mol_rdkit)
        mol["properties"]["isomeric_smiles"] = Chem.MolToSmiles(mol_rdkit, isomericSmiles=True)
        return mol
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
