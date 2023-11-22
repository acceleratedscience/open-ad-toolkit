""" This library is the main library for invoking functions from the command line
it does contain some general functions or plugin calls otherwise invokes toolkits or core / plugin functions"""
#!/usr/local/opt/python@3.9/bin/python3.9
import os
import re
import readline
import webbrowser
import pyperclip
import pandas as pd

# Flask
from openad.flask_apps import launcher
from openad.flask_apps.dataviewer.routes import fetchRoutesDataViewer

import openad.app.login_manager as login_manager

# Core

from openad.core.lang_file_system import import_file, export_file, copy_file, remove_file, list_files
from openad.core.lang_sessions_and_registry import (
    clear_other_sessions,
    write_registry,
    registry_add_toolkit,
    registry_deregister_toolkit,
    initialise_registry,
    update_main_registry_env_var,
)
from openad.core.lang_workspaces import (
    create_workspace,
    remove_workspace,
    list_workspaces,
    set_workspace,
    get_workspace,
)
from openad.core.lang_mols import display_mols
from openad.core.lang_runs import display_run, exec_run, save_run, list_runs
from openad.core.lang_dev import flask_example
from openad.core.grammar import create_statements

# Toolkits
from openad.toolkit.toolkit_main import load_toolkit
from openad.toolkit.toolkit_main import execute_tookit
from openad.toolkit.toolkit_main import load_toolkit_description
from openad.llm_assist.llm_interface import how_do_i, set_llm, clear_llm_auth

# Global variables
from openad.app.global_var_lib import _meta_dir
from openad.app.global_var_lib import _meta_dir_toolkits
from openad.app.global_var_lib import _meta_registry
from openad.app.global_var_lib import _meta_registry_session
from openad.app.global_var_lib import _meta_login_registry
from openad.app.global_var_lib import _meta_workspaces
from openad.app.global_var_lib import _all_toolkits

# Helpers
from openad.helpers.output import (
    msg,
    output_text,
    output_error,
    output_warning,
    output_success,
    output_table,
)
from openad.helpers.general import refresh_prompt, user_input, validate_file_path, ensure_file_path
from openad.helpers.splash import splash
from openad.helpers.output_content import openad_intro

from openad.plugins import edit_json

# Importing our own plugins.
# This is temporary until every plugin is available as a public pypi package.
# from openad.plugins.style_parser import tags_to_markdown


# This is called by the default run_cmd method, for executing current commands.
def lang_parse(cmd_pointer, parser):
    """the routes commands to the correct functions"""
    # Workspace commands
    if parser.getName() == "create_workspace_statement":
        return create_workspace(cmd_pointer, parser)  # Addressed
    elif parser.getName() == "remove_workspace_statement":
        return remove_workspace(cmd_pointer, parser)  # Addressed
    elif parser.getName() == "set_workspace_statement":
        return set_workspace(cmd_pointer, parser)
    elif parser.getName() == "list_workspaces":
        return list_workspaces(cmd_pointer, parser)
    elif parser.getName() == "get_workspace":
        return get_workspace(cmd_pointer, parser)

    # Toolkit commands
    elif parser.getName() == "add_toolkit":
        return registry_add_toolkit(cmd_pointer, parser)
    elif parser.getName() == "remove_toolkit":
        return registry_deregister_toolkit(cmd_pointer, parser)
    elif parser.getName() == "list_toolkits":
        return list_toolkits(cmd_pointer, parser)
    elif parser.getName() == "list_all_toolkits":
        return list_all_toolkits(cmd_pointer, parser)
    elif parser.getName() == "set_context":
        return set_context(cmd_pointer, parser)
    elif parser.getName() == "get_context":
        return get_context(cmd_pointer, parser)
    elif parser.getName() == "unset_context":
        return unset_context(cmd_pointer, parser)
    elif parser.getName() in _all_toolkits:
        # Toolkit welcome screens
        return output_text(splash(parser.getName(), cmd_pointer), nowrap=True)

    # Language Model How To
    elif parser.getName() == "how_do_i":
        result = how_do_i(cmd_pointer, parser)
        if result is False:
            return False
        update_main_registry_env_var(cmd_pointer, "refresh_help_ai", False)
        write_registry(cmd_pointer.settings, cmd_pointer)
        return result
    elif parser.getName() == "set_llm":
        try:
            result = set_llm(cmd_pointer, parser)
            cmd_pointer.llm_model = cmd_pointer.llm_models[cmd_pointer.llm_service]
            update_main_registry_env_var(cmd_pointer, "llm_service", cmd_pointer.llm_service)
            cmd_pointer.refresh_vector = True
            cmd_pointer.refresh_train = True
            write_registry(cmd_pointer.settings, cmd_pointer, False)
            write_registry(cmd_pointer.settings, cmd_pointer, True)
        except (
            Exception  # pylint: disable=broad-exception-caught
        ) as e:  # do not care what exception is, just returning failure
            print(e)
        return result
    elif parser.getName() == "clear_llm_auth":
        result = clear_llm_auth(cmd_pointer, parser)
        return result

    # Run commands
    elif parser.getName() == "create_run":
        # This simply works off the history file entry no procedure needed.
        return output_text(msg("create_run_started"), cmd_pointer, pad=1, nowrap=True)
    elif parser.getName() == "save_run":
        try:
            return save_run(cmd_pointer, parser)
        except (
            Exception  # pylint: disable=broad-exception-caught
        ):  # do not care what exception is, just returning failure
            return False
    elif parser.getName() == "list_runs":
        return list_runs(cmd_pointer, parser)
    elif parser.getName() == "display_run":
        return display_run(cmd_pointer, parser)
    elif parser.getName() == "exec_run":
        exec_run(cmd_pointer, parser)

    # File system commands
    elif parser.getName() == "list_files":
        return list_files(cmd_pointer, parser)
    elif parser.getName() == "import_file":
        return import_file(cmd_pointer, parser)
    elif parser.getName() == "export_file":
        return export_file(cmd_pointer, parser)
    elif parser.getName() == "copy_file":
        return copy_file(cmd_pointer, parser)
    elif parser.getName() == "remove_file":
        return remove_file(cmd_pointer, parser)

    # General commands
    elif parser.getName() == "welcome":
        return output_text(splash(cmd_pointer=cmd_pointer), nowrap=True)
    elif parser.getName() == "get_status":
        return get_status(cmd_pointer, parser)
    elif parser.getName() == "display_history":  # Addressed
        return display_history(cmd_pointer, parser)
    elif parser.getName() == "display_data":
        return display_data(cmd_pointer, parser)
    elif parser.getName() == "display_data__save":
        return display_data__save(cmd_pointer, parser)
    elif parser.getName() == "display_data__open":
        return display_data__open(cmd_pointer, parser)
    elif parser.getName() == "display_data__edit":
        return display_data__open(cmd_pointer, parser, True)
    elif parser.getName() == "display_data__copy":
        return display_data__copy(cmd_pointer, parser)
    elif parser.getName() == "display_data__display":
        return display_data__display(cmd_pointer, parser)
    elif parser.getName() == "clear_sessions":
        return clear_other_sessions(cmd_pointer, parser)
    elif parser.getName() == "edit_config":
        return edit_config(cmd_pointer, parser)

    # Help commands
    elif parser.getName() == "intro":
        return output_text(openad_intro, edge=True, pad=3)
    elif parser.getName() == "docs":
        return docs(cmd_pointer, parser)

    # Show molecules commands
    elif parser.getName() == "show_molecules":
        return display_mols(cmd_pointer, parser)
    elif parser.getName() == "show_molecules_df":
        return display_mols(cmd_pointer, parser)

    # Toolkit execution
    elif str(parser.getName()).startswith("toolkit_exec_"):
        try:
            return execute_tookit(cmd_pointer, parser)
        except (
            Exception  # pylint: disable=broad-exception-caught
        ) as err:  # do not care what exception is, just returning failure
            err = err + "\n" + str(parser.asList())
            return output_error(msg("fail_toolkit_exec_cmd"), cmd_pointer)

    # Development commands (unpublished in help)
    elif parser.getName() == "flask_example":
        return flask_example(cmd_pointer, parser)

    return


# Initialises the metadata and workspace directorys for the tool when first run
# if a directory has been deleted it will recreate it
def initialise():
    """Initialise key paths"""
    if not os.path.isdir(_meta_dir):
        os.mkdir(_meta_dir)
    if not os.path.isdir(_meta_dir_toolkits):
        os.mkdir(_meta_dir_toolkits)
    if not os.path.isdir(_meta_workspaces + "/DEFAULT"):
        os.mkdir(_meta_workspaces)
        os.mkdir(_meta_workspaces + "/DEFAULT")
    if not os.path.isfile(_meta_registry):
        initialise_registry()
    if not os.path.isdir(os.path.dirname(_meta_registry_session)):
        os.mkdir(os.path.dirname(_meta_registry_session))
    if not os.path.isfile(_meta_login_registry):
        login_manager.initialise_toolkit_login()


# Intialize the stateful pickle for the user's sessions.


# Open documentation webpage.
def docs(cmd_pointer, parser):  # pylint: disable=unused-argument # generic pass through used or unused
    """points to online documentation"""
    url = "https://research.ibm.com/topics/accelerated-discovery"
    webbrowser.open_new(url)
    return output_warning(
        "Our documentation website is yet to be built,\nbut this command will open it in the browser.",
        cmd_pointer,
    )


# adds a registry Toolkit Directory
# in future it will take a tar file and explode it into the directory...
# future work on packaging etc...
def welcome(cmd_pointer, parser):  # pylint: disable=unused-argument # generic pass through used or unused
    """Display welcome screen"""
    return output_text(splash(), nowrap=True)


# Display the current context and workspace.
def get_status(cmd_pointer, parser):  # pylint: disable=unused-argument # generic pass through used or unused
    """gets current status and returns to user"""
    status = "\n".join(
        (
            f'<yellow>Current workspace</yellow>: {cmd_pointer.settings["workspace"]}',
            f'<yellow>Current context</yellow>: {str(cmd_pointer.settings["context"])}',
            "<soft>To see more details, run <cmd>get workspace</cmd> or <cmd>get context</cmd>.</soft>",
        )
    )
    return output_text(status, cmd_pointer, nowrap=True, pad=1)


# List the installed toolkits.
def list_toolkits(cmd_pointer, parser):  # pylint: disable=unused-argument # generic pass through used or unused
    """list installed toolkits"""
    toolkits = []
    table_headers = ("Toolkit", "Description")

    # Assemble table data.
    for name in cmd_pointer.settings["toolkits"]:
        description = load_toolkit_description(cmd_pointer, name)
        toolkits.append(list([name, description]))

    # No toolkits installed yet.
    if len(toolkits) == 0:
        return output_text(msg("no_toolkits_installed"), cmd_pointer, pad=1)

    # Display/return table.
    return output_table(toolkits, cmd_pointer, headers=table_headers, note=msg("all_toolkits_currently_installed"))


# List all available toolkits
def list_all_toolkits(cmd_pointer, parser):  # pylint: disable=unused-argument # generic pass through used or unused
    # This will need to be replaced with a scan of the toolkits directory.
    """lists all toolkits"""
    toolkits = []
    table_headers = ("Toolkit", "Installed", "Description")

    # Assemble table data.
    for name in _all_toolkits:
        is_installed = "Yes" if name in cmd_pointer.settings["toolkits"] else "-"
        description = load_toolkit_description(cmd_pointer, name)
        toolkits.append(list([name, is_installed, description]))

    # Add styling tags
    for i, row in enumerate(toolkits):
        is_installed = row[1] == "Yes"
        if not is_installed:
            for j, col_text in enumerate(row):
                toolkits[i][j] = f"<soft>{col_text}</soft>"

    # Display/return table.
    return output_table(toolkits, cmd_pointer, headers=table_headers)


# Set the context of the application to and existing toolkit.
# This means the user will only receive access to base commands
# and the toolkit commands of the toolkit currently in context.
def set_context(cmd_pointer, parser):
    """Sets current toolkit context"""

    reset = False
    if "reset" in parser:
        reset = True

    toolkit_name = parser["toolkit_name"].upper()
    toolkit_current = None

    if toolkit_name.upper() not in cmd_pointer.settings["toolkits"]:
        # if toolkit_name is None: # Trash, without toolkit_name the command is invalidated
        #     return get_context(cmd_pointer, parser)

        # Toolkit doesn't exist.
        return output_error(msg("fail_toolkit_not_installed", toolkit_name, split=True), cmd_pointer, nowrap=True)

    else:
        old_cmd_pointer_context = cmd_pointer.settings["context"]
        old_toolkit_current = cmd_pointer.toolkit_current
        load_ok, toolkit_current = load_toolkit(toolkit_name)

        if load_ok:
            cmd_pointer.settings["context"] = toolkit_name
            cmd_pointer.toolkit_current = toolkit_current
            refresh_prompt(cmd_pointer.settings)
            write_registry(cmd_pointer.settings, cmd_pointer)
            create_statements(cmd_pointer)
            cmd_pointer.current_help.reset_help()
            # cmd_pointer.current_help.help_current.extend(toolkit_current.methods_help)
            login_success = False
            expiry_datetime = None
            try:
                # raise BaseException('Error message') # For testing
                login_success, expiry_datetime = login_manager.load_login_api(cmd_pointer, toolkit_name, reset=reset)

            except (
                Exception  # pylint: disable=broad-exception-caught
            ) as err:  # do not care what exception is, just returning failure
                # Error logging in.
                output_error(msg("err_login", toolkit_name, err, split=True), cmd_pointer=cmd_pointer, return_val=False)
                cmd_pointer.settings["context"] = old_cmd_pointer_context
                cmd_pointer.toolkit_current = old_toolkit_current
                return False

            if not login_success:
                # Failed to log in. # error reporting handled by Toolkit
                cmd_pointer.settings["context"] = old_cmd_pointer_context
                cmd_pointer.toolkit_current = old_toolkit_current
                unset_context(cmd_pointer, None)
                return False

            # Success switching context & loggin in.
            if old_cmd_pointer_context != cmd_pointer.settings["context"]:
                if cmd_pointer.notebook_mode or cmd_pointer.api_mode:
                    return output_success(
                        msg("success_login", toolkit_name, expiry_datetime, split=True),
                        cmd_pointer=cmd_pointer,
                        return_val=False,
                    )
                else:
                    return output_text(splash(toolkit_name, cmd_pointer), nowrap=True)

        else:
            # Failed to load the toolkit
            cmd_pointer.settings["context"] = old_cmd_pointer_context
            cmd_pointer.toolkit_current = old_toolkit_current
            return output_error(msg("fail_toolkit_loading", toolkit_name), cmd_pointer)


# Display your current context.
def get_context(cmd_pointer, parser):  # pylint: disable=unused-argument # generic pass through used or unused
    """gets current toolkit context"""
    current_toolkit = cmd_pointer.settings["context"]
    return output_text(splash(current_toolkit, cmd_pointer), nowrap=True)


# Unset the context of the application.
def unset_context(cmd_pointer, parser):  # pylint: disable=unused-argument # generic pass through used or unused
    """Unsets current toolkit Context"""
    if cmd_pointer.settings["context"] is None:
        return output_text(msg("no_context_set"), cmd_pointer, pad=1)
    cmd_pointer.settings["context"] = None
    cmd_pointer.toolkit_current = None
    cmd_pointer.current_help.reset_help()
    write_registry(cmd_pointer.settings, cmd_pointer)
    refresh_prompt(cmd_pointer.settings)
    create_statements(cmd_pointer)


# Display history of commands.
def display_history(cmd_pointer, parser):  # pylint: disable=unused-argument # generic pass through used or unused
    """Displays last 30 items in current workspace history file"""
    history = []

    # Fetch last 30 commands from history.
    i = 0
    hist_len = readline.get_current_history_length()
    index_col_width = 3
    min_gap = 2
    reached_bottom = False

    while i < 30:
        # Make sure all items are aligned.
        line_nr = hist_len - i
        line_nr_length = len(str(line_nr))
        if line_nr_length + min_gap > index_col_width:
            index_col_width = line_nr_length + min_gap
        # gap = (index_col_width - line_nr_length) * " " no longer requred

        # Fetch history item.
        try:
            entry = readline.get_history_item(line_nr)
            i = i + 1
            if entry:
                entry = entry.replace("\n", "")
                entry_str = f"<cmd>{entry}</cmd>"
            else:
                entry_str = "<soft>Workspace created</soft>"
            history.append((line_nr, entry_str))

            # Reached bottom.
            if not entry:
                reached_bottom = True
                break
        except (
            Exception  # pylint: disable=broad-exception-caught
            # do not care what exception is, just returning failure
        ) as err:
            output_error(msg("err_fetch_history", err, split=True), cmd_pointer)
            i = 31

    # Add ellipsis if history is longer than 30 items.
    if not reached_bottom:
        history.append((line_nr - 1, "<soft>...</soft>"))
    history.reverse()

    # Display/return table.
    return output_table(history, cmd_pointer, headers=["", "Command History"])


# Display a csv file in a table.
def display_data(cmd_pointer, parser):
    """display data in a file"""
    workspace_path = cmd_pointer.workspace_path(cmd_pointer.settings["workspace"].upper()) + "/"
    file_path = parser["file_path"]
    filename = file_path.split("/")[-1]

    # Allow for no extension.
    if len(filename.split(".")) == 1:
        filename = filename + ".csv"
        file_path = file_path + ".csv"

    # Open file
    try:
        if filename.split(".")[-1].lower() == "csv":
            # From csv file.
            try:
                df = pd.read_csv(workspace_path + file_path)
                return output_table(df, cmd_pointer, is_data=True)
            except FileNotFoundError:
                return output_error(msg("fail_file_doesnt_exist", file_path), cmd_pointer)
            except Exception as err:  # pylint: disable=broad-exception-caught
                # do not care what exception is, just returning failure
                return output_error(msg("err_load_csv", err, split=True), cmd_pointer)
        else:
            # Other file formats --> error.
            return output_error(msg("invalid_file_format", "csv", split=True), cmd_pointer)

    except Exception as err:  # pylint: disable=broad-exception-caught
        output_error(msg("err_unknown", err, split=True), cmd_pointer)


# --> Save data to a csv file.
def display_data__save(cmd_pointer, parser):
    """saves data from viewer"""
    # Preserve memory for further follow-up commands.
    cmd_pointer.memory.preserve()

    data = cmd_pointer.memory.get()
    if data is None:
        return output_error(msg("memory_empty", "display"), pad=1)

    # Set variables.
    workspace_path = cmd_pointer.workspace_path(cmd_pointer.settings["workspace"].upper()) + "/"
    file_path = parser["file_path"] if "file_path" in parser else None  # Parser as_dict?
    file_path = validate_file_path(file_path, ["csv"], cmd_pointer)

    # Prompt file path if missing.
    while not file_path:
        file_path = user_input(cmd_pointer, "Filename")
        file_path = validate_file_path(file_path, ["csv"], cmd_pointer)

    # Ensure the file_path is kosher:
    # - Make sure we won't override an existing file
    # - Create folder structure if it doesn't exist yet
    file_path_ok = ensure_file_path(workspace_path + file_path)

    # Save data to file.
    if file_path_ok:
        data.to_csv(workspace_path + file_path, index=False)
        return output_success(msg("success_save_data", file_path), cmd_pointer)
    else:
        return output_error(msg("fail_save_data"), cmd_pointer)


# --> Open data in browser UI.
def display_data__open(
    cmd_pointer, parser, edit_mode=False
):  # pylint: disable=unused-argument # generic pass through used or unused
    # Preserve memory for further follow-up commands.
    """open display data"""
    cmd_pointer.memory.preserve()

    data = cmd_pointer.memory.get()
    if data is None:
        return output_error(msg("memory_empty", "display"), pad=1)

    # Load routes and launch browser UI.
    data = data.to_json(orient="records")
    routes = fetchRoutesDataViewer(data, cmd_pointer)
    a_hash = "#edit" if edit_mode else ""
    launcher.launch(cmd_pointer, routes, "dataviewer", hash=a_hash)


# --> Open data in browser UI.
def display_data__copy(cmd_pointer, parser):  # pylint: disable=unused-argument # generic pass through used or unused
    """displays copy of data in data viewer"""
    # Preserve memory for further follow-up commands.
    cmd_pointer.memory.preserve()

    data = cmd_pointer.memory.get()
    if data is None:
        return output_error(msg("memory_empty", "display"), pad=1)

    pyperclip.copy(data.to_csv(sep="\t"))
    return output_text(msg("data_copied"), pad=1)


# --> Display result in the CLI or Jupyter.
def display_data__display(cmd_pointer, parser):  # pylint: disable=unused-argument # generic pass through used or unused
    """displays last result set in viewer"""
    # Preserve memory for further follow-up commands.
    cmd_pointer.memory.preserve()

    data = cmd_pointer.memory.get()
    if data is None:
        return output_error(msg("memory_empty", "display"), pad=1)

    is_df = isinstance(data, pd.DataFrame)
    if is_df:
        return output_table(data, cmd_pointer, is_data=True)
    else:
        return output_text(data, pad=1)


# Edit a JSON config file.
def edit_config(cmd_pointer, parser):
    """Edits a json document in current workspace"""
    workspace_path = cmd_pointer.workspace_path(cmd_pointer.settings["workspace"].upper()) + "/"

    # Abort in Jupyter.
    if cmd_pointer.notebook_mode is True:
        print("Editing JSON files is only available from command line.")
        return True

    # Load schema.
    schema = None
    if "schema" in parser.as_dict():
        # From parameter.
        schema = workspace_path + parser.as_dict()["schema"]
    else:
        # Scan for same filename with -schema suffix.

        schema = workspace_path + re.sub(r"\.json$", "", parser.as_dict()["json_file"]) + "-schema.json"
        if not os.path.isfile(schema):
            schema = None

    # Load JSON file.
    file_to_edit = workspace_path + parser.as_dict()["json_file"]
    if os.path.isfile(file_to_edit):
        edit_json(file_to_edit, schema)  # pylint: disable=not-callable # it is callable
    elif schema:
        # JSON file not found, create new from schema.
        edit_json(file_to_edit, schema, new=True)  # pylint: disable=not-callable # it is callable
    else:
        return output_error(msg("fail_file_doesnt_exist", parser.as_dict()["json_file"]), cmd_pointer)

    return True
