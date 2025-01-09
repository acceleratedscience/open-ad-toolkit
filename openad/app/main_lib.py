"""This library is the main library for invoking functions from the command line
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
from openad.gui.gui_launcher import gui_init
from openad.flask_apps.dataviewer.routes import fetchRoutesDataViewer
from openad.openad_model_plugin.openad_model_toolkit import openad_model_requestor
from openad.openad_model_plugin.catalog_model_services import (
    catalog_add_model_service,
    uncatalog_model_service,
    get_catalog_namespaces,
    model_service_status,
    service_down,
    service_up,
    local_service_up,
    model_service_config,
    add_service_auth_group,
    remove_service_auth_group,
    attach_service_auth_group,
    detach_service_auth_group,
    list_auth_services,
    get_model_service_result,
)

# molecules
from openad.smols.smol_functions import df_has_molecules
from openad.smols.smol_batch_files import load_mols_to_mws, merge_molecule_property_data
from openad.smols.smol_commands import (
    display_molecule,
    display_property_sources,
    add_molecule,
    remove_molecule,
    list_molecules,
    show_molecules,
    save_molecules_DEPRECATED,
    load_molecules_DEPRECATED,
    display_molsets_DEPRECATED,
    export_molecule,
    get_smol_prop,
    get_smol_prop_lookup_error,
    rename_mol_in_list,
    clear_molecules,
    export_mws,
    show_mol,
    show_molset,
    show_molset_df,
    merge_molecules_DEPRECATED,
)

from openad.mmols.mmol_commands import show_mmol
from openad.smols.smol_cache import enrich_mws_with_analysis, clear_analysis
import openad.app.login_manager as login_manager

# Core
from openad.core.lang_file_system import import_file, export_file, copy_file, remove_file, open_file, list_files
from openad.core.lang_sessions_and_registry import (
    clear_sessions,
    write_registry,
    registry_add_toolkit,
    registry_remove_toolkit,
    update_toolkit,
    update_all_toolkits,
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
from openad.core.lang_runs import display_run, exec_run, save_run, list_runs, remove_run
from openad.core.lang_dev import flask_example
from openad.core.grammar import create_statements

# Toolkits
from openad.toolkit.toolkit_main import load_toolkit
from openad.toolkit.toolkit_main import execute_tookit
from openad.toolkit.toolkit_main import load_toolkit_description
from openad.llm_assist.llm_interface import how_do_i, set_llm, clear_llm_auth

# GUI
from openad.gui.gui_commands import install_gui, launch_gui, restart_gui, quit_gui

# Global variables
from openad.app.global_var_lib import _meta_dir
from openad.app.global_var_lib import _meta_dir_toolkits
from openad.app.global_var_lib import _meta_registry
from openad.app.global_var_lib import _meta_registry_session
from openad.app.global_var_lib import _meta_login_registry
from openad.app.global_var_lib import _meta_workspaces
from openad.app.global_var_lib import _all_toolkits
from openad.app.global_var_lib import MEMORY
from openad.app.global_var_lib import GLOBAL_SETTINGS

# Helpers
from openad.helpers.output import output_text, output_error, output_success, output_table
from openad.helpers.output_msgs import msg
from openad.helpers.general import refresh_prompt, user_input, validate_file_path, ensure_file_path
from openad.helpers.splash import splash
from openad.helpers.output_content import openad_intro
from openad.helpers.plugins import display_plugin_overview

from openad.plugins import edit_json

# Importing our own plugins.
# This is temporary until every plugin is available as a public pypi package.
# from openad.plugins.style_parser import tags_to_markdown


# This is called by the default run_cmd method, for executing current commands.
def lang_parse(cmd_pointer, parser):
    """the routes commands to the correct functions"""

    # print("lang_parse", parser.getName())
    # print(parser)

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
        return registry_remove_toolkit(cmd_pointer, parser)
    elif parser.getName() == "update_toolkit":
        return update_toolkit(cmd_pointer, parser)
    elif parser.getName() == "update_all_toolkits":
        return update_all_toolkits(cmd_pointer, parser)
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

    # Model Service grammar
    elif parser.getName() == "get_model_service_result":
        return get_model_service_result(cmd_pointer, parser)
    elif parser.getName() == "catalog_add_model_service":
        result = catalog_add_model_service(cmd_pointer, parser)
        if result is True:
            # update grammer new service added
            create_statements(cmd_pointer)
        return result
    elif parser.getName() == "uncatalog_model_service":
        result = uncatalog_model_service(cmd_pointer, parser)
        if result is True:
            # update grammer service removed
            create_statements(cmd_pointer)
        return result

    elif parser.getName() == "model_service_status":
        # update grammer for definitions not fetched because service was down
        create_statements(cmd_pointer)
        return model_service_status(cmd_pointer, parser)
    elif parser.getName() == "model_service_config":
        return model_service_config(cmd_pointer, parser)
    elif parser.getName() == "get_catalog_namespaces":
        return get_catalog_namespaces(cmd_pointer, parser)
    elif parser.getName() == "service_up":
        return service_up(cmd_pointer, parser)
    elif parser.getName() == "local_service_up":
        return local_service_up(cmd_pointer, parser)
    elif parser.getName() == "service_down":
        return service_down(cmd_pointer, parser)
    elif parser.getName() == "add_service_auth_group":
        return add_service_auth_group(cmd_pointer, parser)
    elif parser.getName() == "remove_service_auth_group":
        return remove_service_auth_group(cmd_pointer, parser)
    elif parser.getName() == "attach_service_auth_group":
        return attach_service_auth_group(cmd_pointer, parser)
    elif parser.getName() == "detach_service_auth_group":
        return detach_service_auth_group(cmd_pointer, parser)
    elif parser.getName() == "list_auth_services":
        return list_auth_services(cmd_pointer, parser)

    # @later -- move out all logic from here.
    # Language Model How To
    elif parser.getName() == "how_do_i":
        result = how_do_i(cmd_pointer, parser)
        if result is False:
            return False
        cmd_pointer.settings["env_vars"]["refresh_help_ai"] = False
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
            cmd_pointer.settings["env_vars"]["refresh_help_ai"] = True
            write_registry(cmd_pointer.settings, cmd_pointer, False)
            write_registry(cmd_pointer.settings, cmd_pointer, True)
            # update_main_registry_env_var(cmd_pointer, "refresh_help_ai", True)
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
        return output_text(msg("create_run_started"), pad=1, nowrap=True)
    elif parser.getName() == "save_run":
        try:
            return save_run(cmd_pointer, parser)
        except (
            Exception  # pylint: disable=broad-exception-caught
        ):  # do not care what exception is, just returning failure
            return False
    elif parser.getName() == "list_runs":
        return list_runs(cmd_pointer, parser)
    elif parser.getName() == "remove_run":
        return remove_run(cmd_pointer, parser)
    elif parser.getName() == "display_run":
        return display_run(cmd_pointer, parser)
    elif parser.getName() == "exec_run":
        exec_run(cmd_pointer, parser)

    # Molecules
    elif parser.getName() == "display_molecule":
        return display_molecule(cmd_pointer, parser)
    elif parser.getName() == "display_property_sources":
        return display_property_sources(cmd_pointer, parser)
    elif parser.getName() == "add_molecule":
        return add_molecule(cmd_pointer, parser)
    elif parser.getName() == "remove_molecule":
        return remove_molecule(cmd_pointer, parser)
    elif parser.getName() == "list_molecules":
        return list_molecules(cmd_pointer, parser)
    elif parser.getName() == "show_molecules":
        return show_molecules(cmd_pointer, parser)
    elif parser.getName() == "save_molecules_DEPRECATED":
        return save_molecules_DEPRECATED(cmd_pointer, parser)
    elif parser.getName() == "load_molecules_DEPRECATED":
        return load_molecules_DEPRECATED(cmd_pointer, parser)
    elif parser.getName() == "merge_molecules_DEPRECATED":
        return merge_molecules_DEPRECATED(cmd_pointer, parser)
    elif parser.getName() == "list_molecule_sets_DEPRECATED":
        return display_molsets_DEPRECATED(cmd_pointer, parser)
    elif parser.getName() == "enrich_mws_with_analysis":
        return enrich_mws_with_analysis(cmd_pointer, parser)
    elif parser.getName() == "export_molecule":
        return export_molecule(cmd_pointer, parser)
    elif parser.getName() == "clear_analysis":
        return clear_analysis(cmd_pointer, parser)
    elif parser.getName() == "get_smol_prop":
        return get_smol_prop(cmd_pointer, parser)
    elif parser.getName() == "get_smol_prop_lookup_error":
        return get_smol_prop_lookup_error(cmd_pointer, parser)
    elif parser.getName() == "rename_molecule":
        return rename_mol_in_list(cmd_pointer, parser)
    elif parser.getName() == "clear_molecules":
        return clear_molecules(cmd_pointer, parser)
    elif parser.getName() in ["load_molecules_file-DEPRECATED", "load_molecules_dataframe-DEPRECATED"]:
        # MAJOR-RELEASE-TODO
        # Un-comment this in next major release to display deprecation message.
        # output_text(
        #     [
        #         "<on_red> This command is deprecated </on_red>",
        #         "Wrong:   <cmd>load molecules <red>using</red> ... <red>merge with pubchem</red></cmd>",
        #         "Correct: <cmd>load molecules <green>from</green> ... <green>enrich</green></cmd>",
        #     ],
        # )
        return load_mols_to_mws(cmd_pointer, parser)
    elif parser.getName() in ["load_molecules_file", "load_molecules_dataframe"]:
        return load_mols_to_mws(cmd_pointer, parser)
    elif parser.getName() in ["merge_molecules_data_file-DEPRECATED", "merge_molecules_data_dataframe-DEPRECATED"]:
        # MAJOR-RELEASE-TODO
        # Un-comment this in next major release to display deprecation message.
        # output_text(
        #     [
        #         "<on_red> This command is deprecated </on_red>",
        #         "Wrong:   <cmd>merge molecules <red>using</red> ... <red>merge with pubchem</red></cmd>",
        #         "Correct: <cmd>merge molecules <green>from</green> ... <green>enrich</green></cmd>",
        #     ],
        # )
        return merge_molecule_property_data(cmd_pointer, parser)
    elif parser.getName() in ["merge_molecules_data_file", "merge_molecules_data_dataframe"]:
        # NOTE: merge_molecules_data_file is not implemented
        return merge_molecule_property_data(cmd_pointer, parser)
    elif parser.getName() == "export_mws":
        return export_mws(cmd_pointer, parser)
    elif parser.getName() == "show_mol":
        return show_mol(cmd_pointer, parser)
    elif parser.getName() == "show_molset":
        return show_molset(cmd_pointer, parser)
    elif parser.getName() == "show_molset_df":
        return show_molset_df(cmd_pointer, parser)

    # Macromolecules
    elif parser.getName() == "show_mmol":
        return show_mmol(cmd_pointer, parser)

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
    elif parser.getName() == "open_file":
        return open_file(cmd_pointer, parser)

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
    elif parser.getName() == "display_data__as_dataframe":
        return display_data__as_dataframe(cmd_pointer, parser)
    elif parser.getName() == "show_data":
        return show_data(cmd_pointer, parser)
    elif parser.getName() == "clear_sessions":
        return clear_sessions(cmd_pointer, parser)
    elif parser.getName() == "edit_config":
        return edit_config(cmd_pointer, parser)

    # GUI commands
    elif parser.getName() == "install_gui":
        return install_gui(cmd_pointer, parser)
    elif parser.getName() == "launch_gui":
        return launch_gui(cmd_pointer, parser)
    elif parser.getName() == "restart_gui":
        return restart_gui(cmd_pointer, parser)
    elif parser.getName() == "quit_gui":
        return quit_gui(cmd_pointer, parser)

    # Help commands
    elif parser.getName() == "intro":
        return output_text(openad_intro, edge=True, pad=3)
    elif parser.getName() == "docs":
        return docs(cmd_pointer, parser)

    elif "@" in parser.getName() and parser.getName().split("@")[1] in [
        "get_molecule_property",
        "get_crystal_property",
        "get_protein_property",
        "generate_data",
    ]:
        return openad_model_requestor(cmd_pointer, parser)

    # Toolkit execution
    elif str(parser.getName()).startswith("toolkit_exec_"):
        try:
            return execute_tookit(cmd_pointer, parser)
        except (
            Exception  # pylint: disable=broad-exception-caught
        ) as err:  # do not care what exception is, just returning failure
            err = err + "\n" + str(parser.asList())
            return output_error(msg("fail_toolkit_exec_cmd"))

    # Plugin overview screens (name or namspace)
    elif parser.getName().lower() in cmd_pointer.plugins_metadata.keys():
        return display_plugin_overview(cmd_pointer.plugins_metadata[parser.getName().lower()])

    # Toolkit overview screens
    elif parser.getName().upper() in _all_toolkits:
        return output_text(splash(parser.getName(), cmd_pointer), nowrap=True)

    # Plugin commands
    elif parser.getName() in cmd_pointer.plugin_objects.keys():
        return cmd_pointer.plugin_objects[parser.getName()].exec_command(cmd_pointer, parser)

    # Development commands (unpublished in help)
    elif parser.getName() == "flask_example":
        return flask_example(cmd_pointer, parser)
    elif parser.getName() == "cmd_pointer":
        return cmd_pointer

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
    url = "https://acceleratedscience.github.io/openad-docs/commands.html"
    webbrowser.open_new(url)


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
    return output_text(status, nowrap=True, pad=1)


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
        return output_text(msg("no_toolkits_installed"), pad=1)

    # Display/return table.
    return output_table(toolkits, is_data=False, headers=table_headers, note=msg("all_toolkits_currently_installed"))


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
    return output_table(toolkits, is_data=False, headers=table_headers)


# Set the context of the application to and existing toolkit.
# This means the user will only receive access to base commands
# and the toolkit commands of the toolkit currently in context.
def set_context(cmd_pointer, parser):
    """Sets current toolkit context"""

    reset = False
    if "reset" in parser:
        reset = True

    toolkit_name = parser["toolkit_name"].upper()
    set_context_by_name(cmd_pointer, toolkit_name, reset)


# Programatically set the context.
# This is used by the main `set context xyz` command, but also
# by a few other commands like `update context xyz` and `add toolkit xyz`.
def set_context_by_name(cmd_pointer, toolkit_name, reset=False, suppress_splash=False):
    toolkit_current = None

    # Toolkit doesn't exist.
    if toolkit_name.upper() not in cmd_pointer.settings["toolkits"]:
        return output_error(msg("fail_toolkit_not_installed", toolkit_name), nowrap=True)

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

            except Exception as err:  # pylint: disable=broad-exception-caught
                # Error loading login API.
                _handle_login_error(cmd_pointer, toolkit_name, old_toolkit_current, old_cmd_pointer_context, err)
                return False

            if not login_success:
                # Failed to log in.
                err = expiry_datetime  # On fail, error is passed as second parameter instead of expiry.
                _handle_login_error(cmd_pointer, toolkit_name, old_toolkit_current, old_cmd_pointer_context, err)
                return False

            # Success switching context & loggin in.
            if old_cmd_pointer_context != cmd_pointer.settings["context"] and not suppress_splash:
                if GLOBAL_SETTINGS["display"] == "terminal" or GLOBAL_SETTINGS["display"] is None:
                    return output_text(splash(toolkit_name, cmd_pointer), nowrap=True)
                else:
                    return output_success(msg("success_login", toolkit_name, expiry_datetime), return_val=False)
        else:
            # Failed to load the toolkit
            cmd_pointer.settings["context"] = old_cmd_pointer_context
            cmd_pointer.toolkit_current = old_toolkit_current
            return output_error(msg("err_load_toolkit", toolkit_name))


# Handle toolkit login error.
def _handle_login_error(cmd_pointer, toolkit_name, old_toolkit_current, old_cmd_pointer_context, err):
    output_error(msg("err_login", toolkit_name, err), return_val=False)
    cmd_pointer.settings["context"] = old_cmd_pointer_context
    cmd_pointer.toolkit_current = old_toolkit_current
    unset_context(cmd_pointer, None)


# Display your current context.
def get_context(cmd_pointer, parser):  # pylint: disable=unused-argument # generic pass through used or unused
    """gets current toolkit context"""
    current_toolkit = cmd_pointer.settings["context"]
    return output_text(splash(current_toolkit, cmd_pointer), nowrap=True)


# Unset the context of the application.
def unset_context(cmd_pointer, parser):  # pylint: disable=unused-argument # generic pass through used or unused
    """Unsets current toolkit Context"""
    if cmd_pointer.settings["context"] is None:
        return output_text(msg("no_context_set"), pad=1)
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
            output_error(msg("err_fetch_history", err))
            i = 31

    # Add ellipsis if history is longer than 30 items.
    if not reached_bottom:
        history.append((line_nr - 1, "<soft>...</soft>"))
    history.reverse()

    # Display/return table.
    return output_table(history, is_data=False, headers=["", "Command History"])


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
                df = df.fillna("")  # Replace NaN with empty string
                return output_table(df)
            except FileNotFoundError:
                return output_error(msg("err_file_doesnt_exist", file_path))
            except Exception as err:  # pylint: disable=broad-exception-caught
                # do not care what exception is, just returning failure
                return output_error(msg("err_load", "CSV", err))
        else:
            # Other file formats --> error.
            return output_error(msg("err_invalid_file_format", "csv"))

    except Exception as err:  # pylint: disable=broad-exception-caught
        output_error(msg("err_unknown", err))


# --> Save data to a csv file.
def display_data__save(cmd_pointer, parser):
    """saves data from viewer"""
    # Preserve memory for further follow-up commands.
    MEMORY.preserve()

    data = MEMORY.get()
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
        return output_success(msg("success_save_data", file_path))
    else:
        return output_error(msg("err_save_data"))


# --> Open data in browser UI.
def display_data__open(
    cmd_pointer, parser, edit_mode=False
):  # pylint: disable=unused-argument # generic pass through used or unused
    """open display data"""
    # Preserve memory for further follow-up commands.
    MEMORY.preserve()

    df = MEMORY.get()
    if df is None:
        return output_error(msg("memory_empty", "display"), pad=1)

    # If there's molecules in the dataframe, open the result in the molset viewer.

    if df_has_molecules(df) and "as_data" not in parser:
        gui_init(cmd_pointer, "result")
        return

    # Once the dataviewer is integrated, this will all be handled by the results page.
    # but until then, when no molecules are detected we spin up the legacy flask dataviewer.

    # Load routes and launch browser UI.
    df = df.to_json(orient="records")
    routes = fetchRoutesDataViewer(df)
    a_hash = "#edit" if edit_mode else ""
    return launcher.launch(cmd_pointer, routes, "dataviewer", hash=a_hash)


# --> Open data in browser UI.
def display_data__copy(cmd_pointer, parser):  # pylint: disable=unused-argument # generic pass through used or unused
    """displays copy of data in data viewer"""
    # Preserve memory for further follow-up commands.
    MEMORY.preserve()

    data = MEMORY.get()
    if data is None:
        return output_error(msg("memory_empty", "display"), pad=1)

    pyperclip.copy(data.to_csv(sep="\t"))
    return output_text(msg("data_copied"), pad=1)


# --> Display result in the CLI or Jupyter.
def display_data__display(cmd_pointer, parser):  # pylint: disable=unused-argument # generic pass through used or unused
    """displays last result set in viewer"""
    # Preserve memory for further follow-up commands.
    MEMORY.preserve()

    data = MEMORY.get()
    if data is None:
        return output_error(msg("memory_empty", "display"), pad=1)

    is_df = isinstance(data, pd.DataFrame)
    if is_df:
        return output_table(data)
    else:
        return output_text(data, pad=1)


# --> Return result as dataframe
def display_data__as_dataframe(cmd_pointer, parser):  # pylint: disable=unused-argument
    """displays last result set in viewer"""
    # Preserve memory for further follow-up commands.
    MEMORY.preserve()

    data = MEMORY.get()
    if data is None:
        return None

    is_df = isinstance(data, pd.DataFrame)
    if is_df:
        return data
    else:
        return None


def show_data(cmd_pointer, parser):
    """
    Show data in GUI.
    """

    file_path = parser["file_path"]
    filename = file_path.split("/")[-1]

    # Allow for no extension.
    if len(filename.split(".")) == 1:
        file_path = file_path + ".csv"

    gui_init(cmd_pointer, "~/" + file_path)


# Edit a JSON config file.
def edit_config(cmd_pointer, parser):
    """Edits a json document in current workspace"""
    workspace_path = cmd_pointer.workspace_path(cmd_pointer.settings["workspace"].upper()) + "/"

    # Abort in Jupyter.
    if GLOBAL_SETTINGS["display"] == "notebook":
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
        return output_error(msg("err_file_doesnt_exist", parser.as_dict()["json_file"]))

    return True
