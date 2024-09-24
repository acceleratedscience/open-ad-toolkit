# Trash - part of old flask app


def get_molecule_data(cmd_pointer, molecule_identifier):
    """
    Fetch full enriched molecule from RDKit.
    Used to enrich the molecule data from the molecule viewer Flask app.
    """

    from openad.smols.smol_commands import get_smol_from_pubchem, get_smol_from_mws

    mol = get_smol_from_mws(cmd_pointer, molecule_identifier)

    if mol is None:
        mol = get_smol_from_pubchem(molecule_identifier)
    if mol is not None:
        cmd_pointer.last_external_molecule = mol.copy()
        return mol
        # return molformat_v2(mol)

    else:
        return None
