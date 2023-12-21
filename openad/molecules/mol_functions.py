"""Implements molecule management functions"""


import pubchempy as pcy
from rdkit import Chem, rdBase

blocker = rdBase.BlockLogs()


MOL_NAME_INDEX = "name"
MOL_SMILES_INDEX = "smiles"
MOL_INCHI_INDEX = "inchi"
MOL_INCHIKEY_INDEX = "inchikey"
MOL_CID_INDEX = "cid"
MOL_SDF_INDEX = "sdf"
MOL_FORMULA = "formula"
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

        # inchikey = Chem.rdinchi.MolToInchikey(rdkit_mol)

        # mol_weight = Chem.Descriptors.ExactMolWt(rdkit_mol)

    except:
        if new_smiles is None:
            return None
    mol = {
        "name": Name,
        "synonyms": {},
        "properties": {},
        "commments": {},
        "analysis": [],
    }
    for i in MOL_PROPERTIES:
        mol["properties"][i] = None

    mol["properties"]["molecular_weight"] = mol_weight
    mol["properties"]["inchi"] = inchi

    mol["properties"]["inchikey"] = inchikey
    mol["properties"]["canonical_smiles"] = new_smiles

    mol["properties"]["molecular_formula"] = formula
    return mol


def merge_molecule_properties(molecule_dict, mol):
    """merges a molecules property with those from a dictionary"""

    if mol is None:
        return None
    if "ROMol" in molecule_dict:
        del molecule_dict["ROMol"]

    for key in molecule_dict:
        mol["properties"][key] = molecule_dict[key]

    return mol


def valid_smiles(input_molecule) -> bool:
    """determines if a molecule is valid"""
    blocker = rdBase.BlockLogs()
    input_molecule = Chem.MolFromSmiles(input_molecule, sanitize=False)
    if input_molecule is None:
        return False
    else:
        try:
            Chem.SanitizeMol(input_molecule)
        except:
            return False
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
    OPENAD_MOL_DICT = {
        "name": None,
        "synonyms": {},
        "sources": {},
        "properties": {},
        "commments": {},
        "analysis": [],
    }
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
    identifier_dict["multipoles_3d"] = mol["properties"]["multipoles_3d"]
    identifier_dict["pharmacophore_features_3d"] = mol["properties"]["pharmacophore_features_3d"]
    return identifier_dict


if __name__ == "__main__":
    success, molec, cmp = get_mol_from_name("ibuprofen")
    print(success)
    if success:
        print(get_identifiers(molec))
        print(get_properties(molec))
