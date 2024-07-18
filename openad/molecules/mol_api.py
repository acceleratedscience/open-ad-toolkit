# Trash - part of old flask app


def get_molecule_data(cmd_pointer, molecule_identifier):
    """
    Fetch full enriched molecule from RDKit.
    Used to enrich the molecule data from the molecule viewer Flask app.
    """

    from openad.molecules.mol_commands import retrieve_mol, retrieve_mol_from_list

    mol = retrieve_mol_from_list(cmd_pointer, molecule_identifier)

    if mol is None:
        mol = retrieve_mol(molecule_identifier)
    if mol is not None:
        cmd_pointer.last_external_molecule = mol.copy()
        return mol
        # return molformat_v2(mol)

    else:
        return None
