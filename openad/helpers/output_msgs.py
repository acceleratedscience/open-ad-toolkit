"""
This file contains all output messages displayed by OpenAD.

Output messages are stored in four different formats, or a combination of them:
    - String: Returns simple static messages
    - Tuple: Returns messages with multiple lines
    - List: Returns multiple messages**
    - Lambda function: Returns dynamic messages with variables

    ** This is used by error/warning/success messages, eg. to display
       the main error message in red, plus secondary error information in gray.


Usage:

    messages = {
        # Simple string.
        "example_1": "Hello world",

        # Multi-line string.
        "example_2": ("First line", "Second line"),

        # Two separate messages.
        "example_3": ["First message", "Second message"],

        # Dynamic - simple string.
        "example_4": lambda first_name, last_name: f"Hello {first_name} {last_name}",
        
        # Dynamic - multi-line string.
        "example_5": lambda fav_food: ("I love:", fav_food),

        # Dynamic - two separate messages.
        "example_6": lambda err: ["Something went wrong", err],
    }

    msg("example_1")
    --> "Hello world"

    msg("example_2")
    --> "First line\nSecond line"

    msg("example_3")
    --> ["First message", "Second message"]

    msg("example_4", "John", "Doe")
    --> "Hello John Doe"

    msg("example_5", "Pizza")
    --> "I love:\nPizza"

    msg("example_6", "Error 400")
    --> ["Something went wrong", "Error 400"]

"""

# This file uses VSCode regions. Click this line and hit Cmd+(K, 2) to collapse all.
_messages = {
    ##########################################################################
    # region - FILE SYSTEM
    ##########################################################################
    # Success
    "success_import": lambda source_file, workspace_name: f"Imported the file {source_file} to your {workspace_name} workspace",
    "success_export": lambda source_file, workspace_name, dest_file: f"Copied the file {source_file} from your {workspace_name} workspace to {dest_file}",
    "success_copy": lambda source_file, source_workspace_name, dest_workspace_name: f"Copied the file {source_file} from your {source_workspace_name} to your {dest_workspace_name} workspace",
    "success_delete": lambda file_name, workspace_name: f"Deleted the file {file_name} from your {workspace_name} workspace",
    "success_save_data": lambda file_path: f"Your data was successfully stored as <yellow>{file_path}</yellow>",
    "success_file_saved": lambda filename=None: (
        f"Successfully saved <yellow>{filename}</yellow> to your workspace"
        if filename
        else "Successfully saved file to your workspace"
    ),
    "success_file_saved_updated": lambda filename, updated_filename: f"Warning: A file with the name '{filename}' already exists.\n<success>Saved as <yellow>{updated_filename}</yellow> instead.</success>",
    # Warning
    "war_no_filename_provided": lambda default: f"No filename provided, reverting to the default '{default}'",
    # Error
    "err_invalid_file_format": lambda *args: [
        "Invalid file format",
        f'Allowed file types: <yellow>{ "</yellow>, <yellow>".join(args) }</yellow>',
    ],
    "err_invalid_file_format_target": lambda *args: f'You can only save to <yellow>{ "</yellow>, <yellow>".join(args) }</yellow> files',
    "err_file_doesnt_exist": lambda *args: f"The selected file <yellow>{args[0]}</yellow> does not exist",
    "err_path_doesnt_exist": lambda path: f"The path <yellow>{path}</yellow> does not exist",  # ?
    "err_file_not_found": lambda path: [f"File not found", path],
    "err_file_no_permission_read": lambda file_path: [
        f"You don't have permission to open the file <yellow>{file_path.split('/')[-1]}</yellow>",
        file_path,
    ],
    "err_file_no_permission_write": lambda file_path: [
        f"You don't have permission to write to the file <yellow>{file_path.split('/')[-1]}</yellow>",
        file_path,
    ],
    "err_file_is_dir": lambda file_path: [
        f"The file you try to open (<yellow>{file_path.split('/')[-1]}</yellow>) is a directory",
        file_path,
    ],
    "err_decode": lambda file_path: [
        f"Unable to decode the file <yellow>{file_path.split('/')[-1]}</yellow>",
        file_path,
    ],
    "err_io": lambda file_path, err_msg: [
        f"An I/O error occurred while opening the file <yellow>{file_path.split('/')[-1]}</yellow>",
        err_msg,
        file_path,
    ],
    "err_save_data": "No data was stored",
    "err_import": lambda err: ["Import failed", err],
    "err_export": lambda err: ["Export failed", err],
    "err_copy": lambda err: ["Copying file failed", err],
    "err_delete": lambda err: ["Deleting file failed", err],
    "err_load": lambda source_type, err: [f"Unable to load {source_type} file", err],
    # endregion
    ##########################################################################
    # region - WORKSPACES
    ##########################################################################
    # General
    "workspace_description": lambda workspace_name, description, active=False: (
        "<soft>Active workspace:</soft>" if active else None,
        f"<yellow>{workspace_name}</yellow>",
        description,
    ),
    # Negative
    "no_workspace_description": "<soft>This workspace has no description</soft>",
    "no_workspace_files": lambda workspace_name: f"<soft>No files found in the <yellow>{workspace_name}</yellow> workspace</soft>",
    # Success
    "success_workspace_create": lambda workspace_name, error_creating_dir: f"<success>Workspace <yellow>{workspace_name}</yellow> successfully created</success>\n{error_creating_dir}",
    "success_workspace_set": lambda workspace_name: f"You successfully set your workspace to <yellow>{workspace_name}</yellow>",
    "success_workspace_remove": lambda workspace_name: f"Workspace <yellow>{workspace_name}</yellow> successfully removed",
    # Warnings.
    "warn_workspace_already_active": lambda workspace_name: f"You already are in the <yellow>{workspace_name}</yellow> workspace",
    "warn_workspace_folder_already_exists": lambda workspace_name: f"<yellow>Note:</yellow> The workspace directory for {workspace_name} already existed",
    # Error
    "invalid_workpace_destination": lambda workspace_name: f"The destination workspace <yellow>{workspace_name}</yellow> does not exist",
    "invalid_workpace": lambda workspace_name: f"The workspace <yellow>{workspace_name}</yellow> does not exist",
    "fail_workspace_already_exists": lambda workspace_name: f"The workspace <yellow>{workspace_name}</yellow> already exists",
    "fail_workspace_delete_default": "You can't delete the default workspace",
    "err_workspace_description": lambda err=None: ["Something went wrong fetching the workspace description", err],
    "err_workspace_create": lambda err: ["Something went wrong creating the workspace", err],
    # endregion
    ##########################################################################
    # region - RUNS
    ##########################################################################
    # General
    "create_run_started": [
        "<yellow>Recording started</yellow>",
        "Run any number of commands and end with <cmd>save run as <run_name></cmd>",
    ],
    # Negative
    "no_runs_saved_yet": [
        "<yellow>No runs saved yet</yellow>",
        "To create a run, run <cmd>create run</cmd>",
    ],
    # Success
    "success_run_save": "Successfully saved your run",
    # Error
    "fail_run_display": lambda run_name: [
        f'No run named "{run_name}" found',
        "To see available runs, run <cmd>list runs</cmd>",
    ],
    "fail_run_create": "No <cmd>create run</cmd> found in history",
    "fail_run_execute": lambda run_line: f"Unable to execute <cmd>{run_line}</cmd>",
    # endregion
    ##########################################################################
    # region - TOOLKITS
    ##########################################################################
    # General
    "all_toolkits_currently_installed": (
        "These are the toolkits you currently have installed",
        "To see all available toolkits, run <cmd>list all toolkits</cmd>",  # Repeat A1
    ),
    "toolkit_installed": lambda toolkit_name: (
        "<on_green> This toolkit is installed </on_green>",
        f"<soft>To activate this toolkit, run <cmd>set context {toolkit_name.lower()}</cmd></soft>",  # Repeat B1
        f"<soft>To see what you can do, run <cmd>? {toolkit_name.lower()}</cmd></soft>",  # Repeat C1
    ),
    "success_update_toolkit": lambda toolkit_name: f"The <yellow>{toolkit_name}</yellow> toolkit was successfully updated",
    # Negative
    "no_toolkits_installed": [
        "<yellow>No toolkits installed yet</yellow>",
        "To see all available toolkits, run <cmd>list all toolkits</cmd>",  # Repeat A2
    ],
    "toolkit_already_installed": lambda toolkit_name: [
        f"The <reset>{toolkit_name}</reset> toolkit was already installed",
        f"To activate this toolkit, run <cmd>set context {toolkit_name.lower()}</cmd>",  # Repeat B2
        f"To see what you can do, run <cmd>? {toolkit_name.lower()}</cmd>",  # Repeat C2
    ],
    "no_context_set": "<soft>No context was set</soft>",
    # Success
    "success_set_context": lambda toolkit_name: f"You successfully changed the toolkit context to <yellow>{toolkit_name}</yellow>",
    "success_toolkit_install": lambda toolkit_name: [
        f"The <yellow>{toolkit_name}</yellow> toolkit was successfully installed",
        f"To see what you can do, run <cmd>? {toolkit_name.lower()}</cmd>",  # Repeat C2
    ],
    "success_instructions_txt": lambda toolkit_name: f"The <yellow>{toolkit_name}</yellow> toolkit's instructions.txt file was successfully updated",
    "success_toolkit_remove": lambda toolkit_name: f"The <yellow>{toolkit_name}</yellow> toolkit was removed",
    "success_update_all_toolkits": lambda list: ("The following toolkits were successfully updated:", list),
    # Error
    "fail_toolkit_exec_cmd": "Failed to execute toolkit command",
    "fail_toolkit_not_installed": lambda toolkit_name: [
        f"<on_red> The <yellow>{toolkit_name}</yellow> toolkit is not installed </on_red>",
        f"To install it, run <cmd>add toolkit {toolkit_name.lower()}</cmd>",
    ],
    "fail_this_toolkit_not_installed": lambda toolkit_name: [
        f"<on_red> This toolkit is not installed </on_red>",
        f"To install it, run <cmd>add toolkit {toolkit_name.lower()}</cmd>",
    ],
    "invalid_toolkit": lambda toolkit_name: [
        f"There is no toolkit named <yellow>{toolkit_name}</yellow> available",
        "Please check your spelling",
    ],
    "err_load_toolkit": lambda *args: f"There was an error loading the <yellow>{args[0]}</yellow> toolkit",
    # "err_load_toolkit_description": lambda *args: f"There was an error loading 'descriptions.txt' for the <yellow>{args[0]}</yellow> toolkit", # trash
    "err_invalid_description_txt": lambda *args: f"The 'descriptions.txt' for the <yellow>{args[0]}</yellow> toolkit should contain the line:\n<yellow>{args[1]}</yellow>",
    "fail_toolkit_not_registered": lambda toolkit_name: f"The <yellow>{toolkit_name}</yellow> toolkit is not currently registered",
    "err_toolkit_install": lambda toolkit_name, err: [
        f"There was an error installing the <yellow>{toolkit_name}</yellow> toolkit",
        err,
    ],
    "err_toolkit_remove": lambda toolkit_name, err: [
        f"There was an error removing the <yellow>{toolkit_name}</yellow> toolkit",
        err,
    ],
    "err_add_command": lambda command, fwd_expr, err: [
        "Unable to load the function " + command,
        "<yellow>" + fwd_expr + "</yellow>",
        err,
    ],
    "err_update_all_toolkits": lambda list: ("The following toolkits failed to update:", list),
    # endregion
    ##########################################################################
    # region - MOLECULE GRID
    ##########################################################################
    # General
    "m2g_tip": (
        "Tip: To select what parameters to display:",
        "<cmd>result.display(**{ 'subset':['NAME'], 'tooltip':['SMILES'] })</cmd>",
        "",
        "For more options, see: https://mols2grid.readthedocs.io/en/latest/notebooks/customization.html",
    ),
    "flask_launch": lambda port: (
        f"<yellow>Launching the web server</yellow> - Press ctrl+c to abort",
        f"<link>http://127.0.0.1:{port}</link>",
    ),
    # Success
    "success_m2g_save": lambda count, file: f"Succefully saved {count} molecules to {file}",
    "success_m2g_select": lambda count: f"You've selected {count} molecules:",
    # Negative
    "no_m2g_name_column": "No name column identifed in data set",
    # Error
    "fail_m2g_smiles_col_missing": "SMILES column missing from input file",
    "fail_m2g_save_jupyter": 'The "save as" functionality is not supported in Jupyter. Instead you can save your selection directly from the grid interface.',
    "fail_render_mols2grid": "Something went wrong in render_mols2grid()",
    "err_m2g_open": lambda err: ["There was an error opening the molecule viewer", err],
    # endregion
    ##########################################################################
    # region - MOLECULE FUNCTIONS
    ##########################################################################
    # Error
    "input_str_missing": "No input string was provided",
    "invalid_sdf": lambda err: ["The selected SDF file is invalid", err],
    "invalid_csv": lambda err: ["The selected CSV file is invalid", err],
    "identifier_missing": lambda file_path: ["No identifier (Inchi or SMILES) found in the provided file", file_path],
    "err_mol_not_on_pubchem": "This molecule is not available on PubChem",
    # endregion
    ##########################################################################
    # region - MOLECULE VIEWER (section to be merged with molecule functions)
    ##########################################################################
    # Error
    "err_molfile_not_found": lambda filename, file_path: [
        f"The requested molecule file {filename} does not exist.",
        file_path,
    ],
    "err_molfile_invalid": lambda filename, file_path: [
        f"The requested molecule file {filename} does not exist.",
        file_path,
    ],
    # endregion
    ##########################################################################
    # region - SESSIONS
    ##########################################################################
    # General
    "confirm_clear_sessions": (
        "Terminate all other sessions?",
        "<soft>Make sure you don't have any crucial processes running in other windows</soft>",
    ),
    "abort_clear_sessions": ["Action aborted", "Run <cmd>clear sessions</cmd> and try again"],
    # Negative
    "no_sessions_cleared": "No sessions have been cleared",
    # Success
    "success_clear_sessions": "All sessions have been cleared",
    # Error
    "err_clear_sessions": lambda err: ["Something went wrong clearing the other sessions", err],
    # endregion
    ##########################################################################
    # region - LOGIN
    ##########################################################################
    # Success
    "success_login_init": "Login registry initialized",
    "success_login": lambda toolkit_name, expiry_datetime: [
        f"You successfully logged in to <yellow>{toolkit_name}</yellow>",
        (
            f"Your access token does not have an expiration date"
            if expiry_datetime == "No Expiry"
            else f"Your access token expires on {expiry_datetime}"
        ),
    ],
    # Error
    "error_login_init": lambda err: ["Something went wrong while initializing the registry", err],
    "err_login": lambda toolkit_name, err=None: [
        f"Something went wrong logging you in to {toolkit_name}.\n<reset>Please check your credentials and run <cmd>set context {toolkit_name} reset </cmd></reset>",
        err,
    ],
    "err_set_context": "There was a problem setting the context, no context was set.",
    # endregion
    ##########################################################################
    # region - LLM
    ##########################################################################
    # Success
    "success_llm_auth_cleared": "Your LLM authentication file has been cleared",
    # Error
    "error_no_llm_auth_file": lambda path: f"No LLM authentication file found at <yellow>{path}</yellow>",
    "fail_llm_auth_cleared": lambda err: ["Your LLM authentication file could not be cleared", err],
    # endregion
    ##########################################################################
    # region - OTHER
    ##########################################################################
    # General
    "enter_to_skip": "<soft>Press enter to skip</soft>",
    "abort": "Aborted",
    "data_copied": "<success>Data copied to clipboard</success>",
    "run_?": "Run <cmd>?</cmd> to list all command options.",
    "csv_to_clipboard": "CSV data copied to clipboard",
    # Negative
    "table_headers_dont_match_columns": lambda headers, col_count: [
        f"The provided headers ({len(headers)}) don't match the number of columns in the data ({col_count})",
        headers,
    ],
    "no_data_memory": "No data stored in memory",
    "table_is_empty": "No data to display",
    "memory_empty": lambda action: f"There is no result to {action}.",
    # Error
    # 'invalid_cmd': 'Not a valid command',
    "err_invalid_cmd": lambda err="": ["Not a valid command", err],
    "err_no_cmds_matching": lambda inp: [f'No commands containing "<yellow>{inp}</yellow>"'],
    "err_no_cmds_starting": lambda inp: [f'No commands starting with "<yellow>{inp}</yellow>"'],
    "err_unknown": lambda err: ["Unknown error", err],
    "err_fetch_history": lambda err: ["There was an error fetching the history", err],
    "err_load_splash": ["Something went wrong loading the splash page", "<yellow>But the show must go on</yellow>"],
    "err_key_exit_before_init": "Keyboard-initiated exit before OpenAD was initialised",
    "err_routes_required": "Routes are required to launch a Flask server",
    "err_cmd_pointer_required": "output_table(): cmd_pointer is required to enable follow-up commands",
    "err_pyparsing": [
        "Unexpected pyparsing error\n<warning>Please screenshot and report circumstance to OpenAD team</warning>",
        "Restart Notebook kernel or application to proceed",
    ],
    # endregion
}


# Procure a display message from output_msgs.py.
def msg(msg_name, *args, custom_messages=None):
    """
    Fetches and formats an output message from the messages dictionary.

    Parameters
    ----------
    msg_name (str):
        The name of the message to fetch.
    args:
        Any number of variables that are required for the lambda function.
    """

    # This function can be imported by individual
    # toolkits with custom _messages parameter set.
    # See DS4SD for an example.
    if custom_messages:
        output = custom_messages[msg_name]
    else:
        output = _messages[msg_name]

    # Lambda function -> execute to get output.
    if callable(output):
        output = output(*args)

    # Tuple -> remove Nones, turn into multi-line string.
    if isinstance(output, tuple):
        output = [x for x in output if x is not None]
        output = "\n".join(output)

    # List -> remove Nones.
    elif isinstance(output, list):
        output = [x for x in output if x is not None]

    # Output can be string or list of strings.
    return output
