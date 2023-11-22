from rdkit.Chem import PandasTools
import rdkit
import mols2grid
import signal
import os
import pandas
from openad.helpers.output import (
    msg,
    output_text,
    output_error,
    output_warning,
    output_success,
    output_table,
)
from openad.helpers.general import parse_path_tree, confirm_prompt
from openad.helpers.spinner import spinner

mol_name_cache = {}
# Flask
from flask import render_template, send_from_directory, request


def fetchRoutesMols2Grid(cmd_pointer, parser):
    # File and directory references.
    workspace_path = cmd_pointer.workspace_path(cmd_pointer.settings["workspace"].upper()) + "/"
    origin_file = parser.as_dict()["moles_file"] if "moles_file" in parser.as_dict() else None
    results_file = parser["results_file"] if "results_file" in parser else None  # Parser as_dict?

    # Validate input.
    if parser.getName() != "show_molecules_df":
        # Origin file doesn't exist.
        if origin_file and not os.path.exists(workspace_path + origin_file):
            return None, output_error(msg("fail_file_doesnt_exist", origin_file), cmd_pointer)

        # Invalid origin file type.
        if origin_file and len(origin_file.strip()) > 0:
            if origin_file.split(".")[-1].lower() not in ["sdf", "csv"]:
                return None, output_error(msg("invalid_file_format", "sdf", "csv", split=True), cmd_pointer)
        else:
            return None, output_error(msg("invalid_file_format", "sdf", "csv", split=True), cmd_pointer)

        # Invalid destination file type.
        if results_file is not None:
            if results_file and len(results_file.strip()) > 0:
                if results_file and results_file.split(".")[-1].lower() not in ["sdf", "csv"]:
                    return None, output_error(msg("invalid_file_format_target", "sdf", "csv", split=True), cmd_pointer)
            else:
                return None, output_error(msg("invalid_file_format_target", "sdf", "csv", split=True), cmd_pointer)

    # Parameters used to initialize our instance.
    # Used with mols2grid.MolGrid()
    # https://mols2grid.readthedocs.io/en/latest/api/molgrid.html#molgrid-object
    # Currently empty, but could be used to customize our default grid.
    m2g_params_init = {
        # 'size': (150, 100),
    }

    # Render the mol2grid and return:
    # - the_mol2grid: The mols2grid object
    # - mol_frame: The dataframe used to render the grid
    def render_mols2grid():
        try:
            if parser.getName() == "show_molecules_df":
                # From dataframe.
                try:
                    name = parser.getName() + "_" + parser.as_dict()["in_dataframe"]
                    mol_frame = cmd_pointer.api_variables[parser.as_dict()["in_dataframe"]]
                    mol_frame = _normalize_mol_df(mol_frame, cmd_pointer)
                    the_mols2grid = mols2grid.MolGrid(mol_frame, name=name, smiles_col="SMILES", **m2g_params_init)
                except BaseException as err:
                    output_error(msg("err_load_dataframe", err, split=True), cmd_pointer, return_val=False)
                    return False
            elif origin_file.split(".")[-1].lower() == "sdf":
                # From sdf file
                try:
                    name = origin_file.split("/")[-1]
                    SDFFile = workspace_path + origin_file

                    mol_frame = PandasTools.LoadSDF(SDFFile)
                    the_mols2grid = mols2grid.MolGrid(mol_frame, name=name, **m2g_params_init)
                except BaseException as err:
                    output_error(msg("err_load_sdf", err, split=True), cmd_pointer, return_val=False)
                    return False
            elif origin_file.split(".")[-1].lower() == "csv":
                # From csv file.
                try:
                    name = origin_file.split("/")[-1]
                    SDFFile = workspace_path + origin_file
                    mol_frame = pandas.read_csv(SDFFile)
                    mol_frame = _normalize_mol_df(mol_frame, cmd_pointer)
                    the_mols2grid = mols2grid.MolGrid(mol_frame, name=name, **m2g_params_init)
                except BaseException as err:
                    output_error(msg("err_load_csv", err, split=True), cmd_pointer, return_val=False)
                    return False
            else:
                # This shouldn't happen because invalid file types are caught above.
                output_error(msg("invalid_file_format", "sdf", "csv", split=True), cmd_pointer, return_val=False)
                return False

            return the_mols2grid, mol_frame

        except BaseException as err:
            output_error(msg("err_m2g_open", err, split=True), cmd_pointer)

    if results_file:
        # In Jupyter "save as" is not allowed, because we don't
        # have a submit button here.
        if cmd_pointer.notebook_mode:
            return None, output_error(msg("fail_m2g_save_jupyter"), cmd_pointer)

        # If the user has specified a non-existing directory path
        # for the result file, we first need to get permission
        # to create the missing dirs as it could be a mistake.
        # TODO: Replace this with helper function - helpers.general.ensure_file_path()
        else:
            path_tree = parse_path_tree(results_file)
            dir_path = ""
            create_missing_dirs = False
            if len(path_tree) > 0:
                dir_path = workspace_path + "/".join(path_tree)
                if not os.path.exists(dir_path):
                    if confirm_prompt("Directory does not exist. Create it?"):
                        create_missing_dirs = True
                    else:
                        return None, output_error(msg("abort"), cmd_pointer)

            # If the user has specified a file that already
            # exists, we need to get permission to overwrite it.
            if os.path.exists(workspace_path + results_file):
                if not confirm_prompt("Destination file already exists. Overwrite?"):
                    return None, output_error(msg("abort"), cmd_pointer)

    # Create the mols2grid object.
    m2g = render_mols2grid()
    if not m2g or not isinstance(m2g, tuple):
        return None, output_error(msg("fail_render_mols2grid"), cmd_pointer)
    the_mols2grid, mol_frame = m2g

    # Render grid in Jupyter.
    if cmd_pointer.notebook_mode:
        if "object" in parser:
            # Return the mols2grid object.
            output_text(msg("m2g_tip"), cmd_pointer)  # force_print=True
            return None, the_mols2grid
        else:
            # Display the grid.
            m2g_params = _compile_default_m2g_params(mol_frame)

            return None, the_mols2grid.display(**m2g_params)

    # Render grid in Flask.
    else:
        # Create list of available parameters which we
        # then display in the molecule selector UI.
        available_params = mol_frame.columns.tolist()
        if "IMG" in available_params:
            available_params.remove("IMG")
        if "mols2grid-id" in available_params:
            available_params.remove("mols2grid-id")

        # If SMILES column is missing from the input file, we abort.
        if "SMILES" not in available_params:
            return None, output_error(msg("fail_m2g_smiles_col_missing"), cmd_pointer)

        #
        # ROUTES
        #

        # @app.route('/')
        def home():
            # Parse URL arguments.
            args = dict(request.args)
            if len(args) > 0:
                # Create mol2grid display parameters object.
                m2g_params = {k: v.split(",") for k, v in args.items()}
            else:
                # No arguments are passed, set default display parameters.
                m2g_params = _compile_default_m2g_params(mol_frame)

            # Render the grid.
            # We use a copy of the m2g_params because the_mols2grid.display modifies it.
            import copy

            m2g_params_copy = copy.deepcopy(m2g_params)
            m2g_instance = the_mols2grid.display(**m2g_params_copy)
            return render_template(
                "/molsviewer/index.html",
                data=m2g_instance.data,
                available_params=available_params,
                m2g_params=m2g_params,
            )

        # @app.route('/exit', methods=['POST'])
        def exit():
            output_error(msg("abort"), cmd_pointer)
            _kill_server()
            return "ok"

        # Submit
        # @app.route('/submit', methods=['POST'])
        def submit():
            import json

            indexes = json.loads(request.data.decode("utf-8"))
            if indexes is not None:
                filtered_df = the_mols2grid.dataframe.iloc[indexes]
            else:
                filtered_df = the_mols2grid.dataframe

            if "results_file" in parser:
                # Store file

                # Create the directories if thet don't exist.
                if create_missing_dirs:
                    os.makedirs(dir_path)

                # Save as SDF file.
                if results_file.split(".")[-1].lower() == "sdf":
                    # Add the ROMol object to the dataframe.
                    # Either by overriding the ROMol column, or by creating it.
                    # Without this, no sdf file can be created.
                    for index, row in filtered_df.iterrows():
                        if "ROMol" not in filtered_df.columns:
                            filtered_df["ROMol"] = None
                        filtered_df.at[index, "ROMol"] = rdkit.Chem.MolFromSmiles(row["SMILES"])

                    PandasTools.WriteSDF(
                        filtered_df,
                        workspace_path + results_file,
                        properties=list(filtered_df.columns),
                    )

                # Save as CSV file.
                elif results_file.split(".")[-1].lower() == "csv":
                    # We remove the RoMol object as it's useless in the CSV file.
                    # It can be regenerated from the SMILES string.
                    filtered_df = filtered_df.drop(["ROMol"], axis=1, errors="ignore")

                    # We remove the img column as it's useless in the CSV file.
                    filtered_df = filtered_df.drop(["IMG"], axis=1, errors="ignore")
                    filtered_df.to_csv(workspace_path + results_file, index=False)

                # Success message
                spinner.stop()
                output_success(msg("success_m2g_save", len(indexes), parser["results_file"]), cmd_pointer)
            else:
                # Display results
                note = None  # 'To see what you can do next, run <cmd>next ?</cmd>'
                spinner.stop()
                output_success(msg("success_m2g_select", len(indexes)), cmd_pointer, pad_top=1)
                output_table(filtered_df, cmd_pointer, note=note)

            # Shut down server
            _kill_server()
            return "ok"

        # NOTE: I tried having a spinning animation here, to indicate that the server is running,
        # but this has proven a bit more challenging that I expected and I have to move on for now.
        # Moenen - june 6 '23

        # ATTEMPT #1 - Using our default spinner, but it doesn't exit elegantly on ctrl+c.
        # - -
        # spinner.start('<link>http://127.0.0.1:5000</link>', no_format=True)

        # ATTEMPT #2 - Home made spinner that exists more elegantlty on ctrl+c,
        # but it blocks the process. Need to look deeper into threading.
        # - -
        # Display loading animation.
        # from helpers.general import loader
        # text = output_text('<link>http://127.0.0.1:5000</link>', cmd_pointer, return_val=True)
        # exit_msg = return output_error(msg('abort'), cmd_pointer, return_val=True)
        # loader(text, ',.-*°*-.', exit_msg=exit_msg, no_format=True, on_abort=_kill_server)

        routes = {
            "/": {"func": home, "method": "GET"},
            "/exit": {"func": exit, "method": "GET"},
            "/submit": {"func": submit, "method": "POST"},
        }

        return routes, the_mols2grid


def _normalize_mol_df(mol_df: pandas.DataFrame, cmd_pointer):
    has_name = False
    contains_name = None

    for i in mol_df.columns:
        # Find the name column.
        if str(i.upper()) == "NAME":
            has_name = True
        if contains_name is None and "NAME" in str(i.upper()):
            contains_name = i

        # Normalize any columns we'll be referring to later.
        if str(i.upper()) == "SMILES":
            mol_df.rename(columns={i: "SMILES"}, inplace=True)
        if str(i.upper()) == "ROMOL":
            mol_df.rename(columns={i: "ROMol"}, inplace=True)
        if str(i.upper()) == "IMG":
            mol_df.rename(columns={i: "IMG"}, inplace=True)

    # Normalize name column.
    if has_name == False and contains_name is not None:
        mol_df.rename(columns={contains_name: "NAME"}, inplace=True)

    # Add names when missing.
    try:
        if "NAME" not in mol_df.columns:
            output_warning(msg("no_m2g_name_column"), cmd_pointer)
            spinner.start("Downloading names")
            mol_df["NAME"] = "unknown"
            for i in mol_df.itertuples():
                # print(_smiles_to_iupac('CCC(COC(=O)CS)(C(=O)C(=O)CS)C(=O)C(=O)CS'))
                mol_df.loc[i.Index, "NAME"] = _smiles_to_iupac(mol_df.loc[i.Index, "SMILES"])
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
    import pubchempy

    if smiles in mol_name_cache:
        return mol_name_cache[smiles]
    try:
        compounds = pubchempy.get_compounds(smiles, namespace="smiles")
        match = compounds[0]
        mol_name_cache[smiles] = str(match)
    except BaseException:
        match = "unknown"
    return str(match)


def _compile_default_m2g_params(mol_frame):
    available_params = mol_frame.columns.tolist()
    if "NAME" in available_params:
        m2g_params = {
            "subset": ["NAME"],
            "tooltip": [x for i, x in enumerate(available_params) if (x not in ["NAME", "mols2grid-id", "IMG"])],
        }
        return m2g_params


def _kill_server():
    os.kill(os.getpid(), signal.SIGINT)
