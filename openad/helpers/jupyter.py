import os
from rdkit import Chem
from IPython.display import display, HTML
from openad.smols.smol_functions import get_mol_rdkit, valid_identifier
from openad.helpers.output_msgs import msg
from openad.helpers.output import output_error, output_warning, output_success, output_text


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

    mol_rdkit = get_mol_rdkit(identifier, identifier_type)
    mol_drawer = Chem.Draw.MolDraw2DSVG(300, 300)  # pylint: disable=no-member
    mol_drawer.DrawMolecule(mol_rdkit)
    mol_drawer.FinishDrawing()
    mol_svg = mol_drawer.GetDrawingText()
    img_html = f'<div style="width:300px; height: 300px; margin: 30px 0; border: solid 1px #ddd; display: inline-block; padding: 32px; position: relative"><div style="position: absolute; top: 8px; left: 8px; font-size: 12px; line-height: 12px; color: #999;">INPUT MOLECULE</div>{mol_svg}</div>'
    display(HTML(img_html))


def save_df_as_csv(cmd_pointer, df, dest_file_path):
    """
    Save a pandas dataframe as a CSV file.

    Parameters
    ----------
    cmd_pointer: Cmd
        The command pointer.
    df: pandas.DataFrame
        The dataframe to save.
    dest_file_path: str
        The destination file path, with your workspace as root
        and an optional .csv extension and leading slash.
        All valid:
        - filename
        - folder1/folder2/filename
        - folder1/folder2/filename.csv
        - /folder1/folder2/filename.csv
    """
    # Remove leading slash
    if dest_file_path.startswith("/"):
        dest_file_path = dest_file_path[1:]

    # Remove any number of ../ from the path to avoid storing
    # files outside the workspace (could be abused)
    while dest_file_path.startswith("../"):
        dest_file_path = dest_file_path.replace("../", "")

    # Ensure CSV extension
    if not dest_file_path.endswith(".csv"):
        dest_file_path = dest_file_path + ".csv"

    # Create destination file path directories if they don't exist
    dirs = dest_file_path.split("/")[:-1]
    workspace_path = cmd_pointer.workspace_path()
    for d in dirs:
        workspace_path += "/" + d
        if not os.path.exists(workspace_path):
            os.makedirs(workspace_path)

    # Find next available filename if the file already exists
    absolute_dest_file_path = cmd_pointer.workspace_path() + "/" + dest_file_path
    base, extension = os.path.splitext(dest_file_path)
    counter = 1
    updated_dest_file_path = None
    while os.path.exists(absolute_dest_file_path):
        updated_dest_file_path = f"{base}-{counter}{extension}"
        absolute_dest_file_path = cmd_pointer.workspace_path() + "/" + updated_dest_file_path
        counter += 1

    # Save the file
    df = df.fillna("")  # Replace NaN with empty string
    df.to_csv(absolute_dest_file_path, index=False)

    # Display success message
    if updated_dest_file_path:
        output_warning(msg("success_file_saved_updated", dest_file_path, updated_dest_file_path), return_val=False)
    else:
        output_success(msg("success_file_saved", dest_file_path), return_val=False)

    # Display hint on how to open it
    output_text(f"<soft>To open it, run <cmd>open '{updated_dest_file_path or dest_file_path}'</cmd></soft>", pad_btm=1)
