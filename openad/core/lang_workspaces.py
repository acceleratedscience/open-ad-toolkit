"""Workspace related Functions"""

import os
from time import sleep

import readline

# Core
from openad.core.lang_sessions_and_registry import write_registry, update_main_registry_env_var

# Global variables
from openad.app.global_var_lib import GLOBAL_SETTINGS

# Helpers
from openad.helpers.output import output_text, output_error, output_warning, output_success, output_table
from openad.helpers.output_msgs import msg
from openad.helpers.general import other_sessions_exist, user_input
from openad.helpers.spinner import spinner


# Sets the current workspace from the fgiven workspaces available
def set_workspace(cmd_pointer, parser):
    """Sets the current Workspace"""
    readline.write_history_file(cmd_pointer.histfile)
    current_workspace_name = cmd_pointer.settings["workspace"].upper()
    new_workspace_name = parser["Workspace_Name"].upper()
    if new_workspace_name not in cmd_pointer.settings["workspaces"]:
        return output_error(msg("invalid_workpace", new_workspace_name))

    elif new_workspace_name == current_workspace_name:
        return output_warning(msg("warn_workspace_already_active", new_workspace_name))
    else:
        cmd_pointer.settings["workspace"] = new_workspace_name
        write_registry(cmd_pointer.settings, cmd_pointer)
        cmd_pointer.histfile = os.path.expanduser(cmd_pointer.workspace_path(new_workspace_name) + "/.cmd_history")
        readline.clear_history()

        try:  # Open history file if not corrupt
            if readline and os.path.exists(cmd_pointer.histfile):
                readline.read_history_file(cmd_pointer.histfile)
        except Exception:
            readline.write_history_file(cmd_pointer.histfile)

        readline.write_history_file(cmd_pointer.histfile)
        return output_success(msg("success_workspace_set", new_workspace_name))


# list the available workspaces....
def list_workspaces(cmd_pointer, parser):
    """Lists all Workspaces"""
    workspaces = []
    table_headers = ("Workspace", "Description")
    note = "To see what you can do with a workspace, run <cmd>? workspace</cmd>."
    current_workspace_name = cmd_pointer.settings["workspace"].upper()

    for name in cmd_pointer.settings["workspaces"]:
        # Format current workspace name.
        if name == current_workspace_name:
            if GLOBAL_SETTINGS["display"] == "notebook":
                name_formatted = f"* {name}"
            else:
                name_formatted = output_text(f"* <green>{name}</green>", return_val=True)
        else:
            name_formatted = name

        # Add 'No description' if no description is available.
        if not cmd_pointer.settings["descriptions"][name]:
            cmd_pointer.settings["descriptions"][name] = output_text(msg("no_workspace_description"), return_val=True)

        workspaces.append(list([name_formatted, cmd_pointer.settings["descriptions"][name]]))

    # Display/return table.
    return output_table(workspaces, is_data=False, headers=table_headers, note=note)


# get the details of a workspace
# needs to be fixed up as workspace metadata plan is built out
def get_workspace(cmd_pointer, parser):
    """gets a workspaces details"""
    if "Workspace_Name" in parser.as_dict():
        workspace_name = parser.as_dict()["Workspace_Name"].upper()
        active = False
    else:
        workspace_name = cmd_pointer.settings["workspace"].upper()
        active = True
    if workspace_name in cmd_pointer.settings["descriptions"]:
        description = cmd_pointer.settings["descriptions"][workspace_name]
    else:
        description = None
    description = description if description else msg("no_workspace_description")

    if workspace_name not in cmd_pointer.settings["workspaces"]:
        return output_error(msg("invalid_workpace", workspace_name))
    else:
        return output_text(msg("workspace_description", workspace_name, description, active), pad=1, edge=True)


# Remove workspace and all its metadata files.
def remove_workspace(cmd_pointer, parser):
    """Removes a registered Workspace from Registry"""

    other_sesh = other_sessions_exist(cmd_pointer)
    if other_sesh is True:
        return

    cmd_pointer.refresh_vector = True
    cmd_pointer.refresh_train = True
    cmd_pointer.settings["env_vars"]["refresh_help_ai"] = True
    update_main_registry_env_var(cmd_pointer, "refresh_help_ai", True)

    workspace_name = parser.as_dict()["Workspace_Name"].upper()
    if workspace_name == "DEFAULT":
        return output_error(msg("fail_workspace_delete_default"))
    if workspace_name in cmd_pointer.settings["workspaces"]:
        cmd_pointer.settings["workspaces"].remove(workspace_name)
        if (
            workspace_name in cmd_pointer.settings["paths"]
        ):  # <-- @Phil added to avoid error, but probably shouldn't be empty.
            cmd_pointer.settings["paths"].pop(workspace_name)
        cmd_pointer.settings["workspace"] = "DEFAULT"
        cmd_pointer.histfile = os.path.expanduser(
            cmd_pointer.workspace_path(cmd_pointer.settings["workspace"].upper()) + "/.cmd_history"
        )
        cmd_pointer.settings["descriptions"].pop(workspace_name)
        write_registry(cmd_pointer.settings, cmd_pointer, True)
        write_registry(cmd_pointer.settings, cmd_pointer)
        return output_success(msg("success_workspace_remove", workspace_name))
    else:
        return output_error(msg("invalid_workpace", workspace_name))


def create_workspace(cmd_pointer, parser):
    """Creates a Workspace"""
    # Make sure existing workspace history file is saved.
    readline.write_history_file(cmd_pointer.histfile)
    cmd_pointer.refresh_vector = True
    cmd_pointer.refresh_train = True

    cmd_pointer.settings["env_vars"]["refresh_help_ai"] = True
    update_main_registry_env_var(cmd_pointer, "refresh_help_ai", True)

    # Abort if other sessions are running.
    other_sesh = other_sessions_exist(cmd_pointer)
    if other_sesh is True:
        return

    # Fetch workspace name.
    workspace_name = parser.as_dict()["Workspace_Name"].upper()

    # Abort if workspace already exists.
    if workspace_name in cmd_pointer.settings["workspaces"]:
        if GLOBAL_SETTINGS["display"] != "api":
            return output_error(msg("fail_workspace_already_exists", workspace_name))

    # Store workspace description.
    try:
        if "proj_desc" in parser.as_dict():
            # From parser.
            description = parser.as_dict()["proj_desc"]
        else:
            # From input.
            output_text(msg("enter_to_skip"), pad_top=1)  # force_print=True
            description = user_input(cmd_pointer, "Workspace description")
            if description is None or len(description.strip()) == 0:
                description = "No workspace description available"

        cmd_pointer.settings["descriptions"][workspace_name] = description
        write_registry(cmd_pointer.settings, cmd_pointer, True)  # Create registry
        write_registry(cmd_pointer.settings, cmd_pointer)  # Create session registry
    except Exception as err:
        return output_error(msg("err_workspace_description", err))

    # Create workspace.
    if "w_path" in parser.as_dict():
        path = parser.as_dict()["w_path"]

        # Expand user path: ~/ --> ../
        # from pathlib import PosixPath
        # path = PosixPath(path).expanduser().resolve() # %%
        path = os.path.expanduser(path)

        if not os.path.exists(path):
            return output_error(msg("err_path_doesnt_exist", path))
        cmd_pointer.settings["paths"][workspace_name] = path
    spinner.start("Creating workspace")
    sleep(0.5)  # Ensure the spinner is displayed for at least a moment.
    cmd_pointer.settings["workspaces"].append(workspace_name)
    cmd_pointer.settings["workspace"] = workspace_name
    cmd_pointer.histfile = os.path.expanduser(
        cmd_pointer.workspace_path(cmd_pointer.settings["workspace"].upper()) + "/.cmd_history"
    )

    error_creating_dir = False
    error_other = False
    try:
        workspace_name = cmd_pointer.settings["workspace"].upper()
        dir_path = os.path.expanduser(cmd_pointer.workspace_path(workspace_name))
        if not os.path.exists(dir_path):
            os.mkdir(dir_path)
        else:
            # This currently happens when you remove a workspace and then try to recreate it.
            # @Phil - we probably should move or archive the workspace folder when removing the workspace.
            os.chdir(dir_path)
            error_creating_dir = msg("warn_workspace_folder_already_exists", workspace_name)
        # Main and session registry writes
        write_registry(cmd_pointer.settings, cmd_pointer, True)
        write_registry(cmd_pointer.settings, cmd_pointer)

        readline.clear_history()
        readline.write_history_file(cmd_pointer.histfile)
        # raise ValueError('This is a test error.\n') @later this causes the app to break permamenently.
    except Exception as err:
        error_other = msg("err_workspace_create", err)

    # Show success/errror message.
    if GLOBAL_SETTINGS["display"] != "api":
        if error_other:
            spinner.fail(output_error(error_other, return_val=True, jup_return_format="plain"))
            spinner.start()
            spinner.stop()
            return
        else:
            add_line = bool(error_creating_dir)
            error_creating_dir = error_creating_dir + "\n" if error_creating_dir else ""
            spinner.succeed(
                output_text(
                    msg("success_workspace_create", workspace_name, error_creating_dir),
                    return_val=True,
                    jup_return_format="plain",
                )
            )
            spinner.start()
            spinner.stop()
            if add_line:
                print("")
            return
