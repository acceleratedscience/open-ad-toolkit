from openad.molecules.mol_commands import retrieve_mol_from_list, retrieve_mol, get_identifiers, get_properties


def get_molecule_data(molecule_identifier, cmd_pointer):
    """Gets displayable Data for a given molecule"""
    return_dict = {}

    mol = retrieve_mol_from_list(cmd_pointer, molecule_identifier)

    if mol is None:
        mol = retrieve_mol(molecule_identifier)
    if mol is not None:
        # cmd_pointer.last_external_molecule = mol.copy()
        return_dict["smiles"] = get_identifiers(mol)["canonical_smiles"]
        return_dict["identifiers"] = get_identifiers(mol)
        return_dict["properties"] = get_properties(mol)
        return_dict["synonyms"] = mol["synonyms"]["Synonym"]
        return_dict["analysis"] = mol["analysis"]
        return return_dict
    else:
        return None
