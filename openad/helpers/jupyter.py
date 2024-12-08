from rdkit import Chem
from IPython.display import display, HTML
from openad.smols.smol_functions import get_mol_rdkit, valid_identifier
from openad.helpers.output import output_error


def jup_display_input_molecule(identifier, identifier_type=None):
    """
    Display the input molecule of a query in the Jupyter Notebook.

    Parameters:
    -----------
    identifier: str
        The molecule identifier
    identifier_type: str
        The molecule identifier type:
        SMILES, InChI, InChIKey
    """

    if not valid_identifier(identifier):
        return

    print(1, identifier, identifier_type)
    mol_rdkit = get_mol_rdkit(identifier, identifier_type)
    print(2, mol_rdkit)
    mol_drawer = Chem.Draw.MolDraw2DSVG(300, 300)  # pylint: disable=no-member
    print(3, mol_drawer)
    mol_drawer.DrawMolecule(mol_rdkit)
    print(4)
    mol_drawer.FinishDrawing()
    print(5)
    mol_svg = mol_drawer.GetDrawingText()
    print(6)
    img_html = f'<div style="width:300px; height: 300px; margin: 30px 0; border: solid 1px #ddd; display: inline-block; padding: 32px; position: relative"><div style="position: absolute; top: 8px; left: 8px; font-size: 12px; line-height: 12px; color: #999;">INPUT MOLECULE</div>{mol_svg}</div>'
    print(7)
    display(HTML(img_html))
