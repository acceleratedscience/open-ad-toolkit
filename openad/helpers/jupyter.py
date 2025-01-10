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
    if mol_rdkit:
        try:
            mol_drawer = Chem.Draw.MolDraw2DSVG(300, 300)  # pylint: disable=no-member
            mol_drawer.DrawMolecule(mol_rdkit)
            mol_drawer.FinishDrawing()
            mol_svg = mol_drawer.GetDrawingText()
            img_html = f'<div style="width:300px; height: 300px; margin: 30px 0; border: solid 1px #ddd; display: inline-block; padding: 32px; position: relative"><div style="position: absolute; top: 8px; left: 8px; font-size: 12px; line-height: 12px; color: #999;">INPUT MOLECULE</div>{mol_svg}</div>'
            display(HTML(img_html))
            # raise Exception("This is a test error.")
        except Exception as err:  # pylint: disable=broad-except
            w_id_type = " with identifier type '" + identifier_type + "'" if identifier_type else ""
            ouput_msg = f"Error in jup_display_input_molecule():\nSomething went wrong displaying the molecule '{identifier}'{w_id_type}."
            output_error([ouput_msg, err], return_val=False)


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


def parse_using_clause(params: list, allowed: list):
    """
    Parse the content of the USING clause from a list of tuples into a dictionary.

    - Input: [["foo", 123], ["bar", "baz"]]
    - Output: {"foo": 123, "bar": "baz"}

    Parameters
    ----------
    params : list
        List of tuples
    allowed : list
        List of allowed keys
    """

    params_dict = {}
    invalid_keys = []

    if not params:
        return params_dict

    for [key, val] in params:
        if key in allowed:
            params_dict[key] = val
        else:
            invalid_keys.append(key)

    if invalid_keys:
        output_warning(
            "Warning: Ignored invalid USING parameters:\n- <error>"
            + ("</error>\n- <error>".join(invalid_keys))
            + "</error>",
            return_val=False,
            pad=1,
        )

    return params_dict
