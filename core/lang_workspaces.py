import os
from time import sleep

import readline
# Core
from core.lang_sessions_and_registry import write_registry,update_main_registry_env_var

# Global variables
from global_var_lib import _meta_dir as _meta_dir
from global_var_lib import _meta_dir_toolkits as _meta_dir_toolkits
from global_var_lib import _meta_registry as _meta_registry
from global_var_lib import _meta_registry_session as _meta_registry_session
from global_var_lib import _meta_login_registry as _meta_login_registry
from global_var_lib import _meta_workspaces as _meta_workspaces
from global_var_lib import _meta_registry_settings as _meta_registry_settings

# Helpers
from helpers.output import msg, output_text, output_error, output_warning, output_success, output_table
from helpers.general import other_sessions_exist, user_input, is_notebook_mode
from helpers.spinner import spinner


# Sets the current workspace from the fgiven workspaces available
def set_workspace(cmd_pointer, parser):
  
    readline.write_history_file(cmd_pointer.histfile)
    current_workspace_name = cmd_pointer.settings['workspace'].upper()
    new_workspace_name = parser['Workspace_Name'].upper()
    if new_workspace_name not in cmd_pointer.settings['workspaces']:
        return output_error(msg('invalid_workpace', new_workspace_name), cmd_pointer)

    elif new_workspace_name == current_workspace_name:
        return output_warning(msg('warn_workspace_already_active', new_workspace_name), cmd_pointer)
    else:
        cmd_pointer.settings['workspace'] = new_workspace_name
        write_registry(cmd_pointer.settings, cmd_pointer)
        cmd_pointer.histfile = os.path.expanduser(cmd_pointer.workspace_path(new_workspace_name) + '/.cmd_history')
        readline.clear_history()

        try:  # Open history file if not corrupt
            if readline and os.path.exists(cmd_pointer.histfile):
                readline.read_history_file(cmd_pointer.histfile)
        except BaseException:
            readline.write_history_file(cmd_pointer.histfile)

        readline.write_history_file(cmd_pointer.histfile)
        return output_success(msg('success_workspace_set', new_workspace_name), cmd_pointer)


# list the available workspaces....
def list_workspaces(cmd_pointer, parser):
    workspaces = []
    table_headers = ('Workspace', 'Description')
    note = 'To see what you can do with a workspace, run <cmd>workspace ?</cmd>.'
    current_workspace_name = cmd_pointer.settings['workspace'].upper()

    for name in cmd_pointer.settings['workspaces']:
        # Format current workspace name.
        if name == current_workspace_name:
            if cmd_pointer.notebook_mode:
                name_formatted = f'* {name}'
            else:
                name_formatted = output_text(f'* <green>{name}</green>', cmd_pointer, return_val=True)
        else:
            name_formatted = name

        # Add 'No description' if no description is available.
        if not cmd_pointer.settings['descriptions'][name]:
            cmd_pointer.settings['descriptions'][name] = output_text(msg('no_workspace_description'), cmd_pointer, return_val=True)

        workspaces.append(list([name_formatted, cmd_pointer.settings['descriptions'][name]]))

    # Display/return table.
    return output_table(workspaces, cmd_pointer, headers=table_headers, note=note)


# get the details of a workspace
# needs to be fixed up as workspace metadata plan is built out
def get_workspace(cmd_pointer, parser):
    workspace_name = parser.as_dict()['Workspace_Name'].upper()
    if workspace_name in cmd_pointer.settings['descriptions']:
        description = cmd_pointer.settings['descriptions'][workspace_name]
    else:
        description = None
    description = description if description else msg('no_workspace_description')

    if workspace_name not in cmd_pointer.settings['workspaces']:
        return output_error(msg('invalid_workpace', workspace_name), cmd_pointer)
    else:
        return output_text(msg('workspace_description', workspace_name, description), cmd_pointer, pad=1, edge=True)


# Remove workspace and all its metadata files.
def remove_workspace(cmd_pointer, parser):
    other, error_msg = other_sessions_exist(cmd_pointer)
    cmd_pointer.refresh_vector=True
    cmd_pointer.refresh_train=True
    cmd_pointer.settings['env_vars']['refresh_help_ai']=True
    update_main_registry_env_var(cmd_pointer,'refresh_help_ai',True)
    if other == True:
        return error_msg

    workspace_name = parser.as_dict()['Workspace_Name'].upper()
    if workspace_name == 'DEFAULT':
        return output_error(msg('fail_workspace_delete_default'), cmd_pointer)
    if workspace_name in cmd_pointer.settings['workspaces']:
        cmd_pointer.settings['workspaces'].remove(workspace_name)
        if workspace_name in cmd_pointer.settings['paths']:  # <-- @Phil added to avoid error, but probably shouldn't be empty.
            cmd_pointer.settings['paths'].pop(workspace_name)
        cmd_pointer.settings['workspace'] = 'DEFAULT'
        cmd_pointer.histfile = os.path.expanduser(cmd_pointer.workspace_path(cmd_pointer.settings['workspace'].upper()) + '/.cmd_history')
        cmd_pointer.settings['descriptions'].pop(workspace_name)
        write_registry(cmd_pointer.settings, cmd_pointer, True)
        write_registry(cmd_pointer.settings, cmd_pointer)
        return output_success(msg('success_workspace_remove', workspace_name), cmd_pointer)
    else:
        return output_error(msg('invalid_workpace', workspace_name), cmd_pointer)


def create_workspace(cmd_pointer, parser):
    # Make sure existing workspace History file is saved.
   
    readline.write_history_file(cmd_pointer.histfile)
    cmd_pointer.refresh_vector=True
    cmd_pointer.refresh_train=True
    
    cmd_pointer.settings['env_vars']['refresh_help_ai']=True
    update_main_registry_env_var(cmd_pointer,'refresh_help_ai',True)
    # Abort if other sessions are running.
    other, error_msg = other_sessions_exist(cmd_pointer)
    if other == True:
        return error_msg

    # Fetch workspace name.
    workspace_name = parser.as_dict()['Workspace_Name'].upper()

    # Abort if workspace already exists.
    if workspace_name in cmd_pointer.settings['workspaces']:
        if cmd_pointer.api_mode == False:
            return output_error(msg('fail_workspace_already_exists', workspace_name), cmd_pointer)

    # Store workspace description.
    try:
        if 'proj_desc' in parser.as_dict():
            # From parser.
            description = parser.as_dict()['proj_desc']
        else:
            # From input.
            output_text(msg('enter_to_skip'), cmd_pointer, pad_top=1)  # force_print=True
            description = user_input(cmd_pointer, 'Workspace description')
        cmd_pointer.settings['descriptions'][workspace_name] = description
        write_registry(cmd_pointer.settings, cmd_pointer, True)  # Create registry
        write_registry(cmd_pointer.settings, cmd_pointer)  # Create session registry
    except BaseException as err:
        return output_error(msg('err_workspace_description', err, split=True), cmd_pointer)

    # Create workspace.
    if 'w_path' in parser.as_dict():
        path = parser.as_dict()['w_path']

        # Expand user path: ~/ --> ../
        # from pathlib import PosixPath
        # path = PosixPath(path).expanduser().resolve() #%%
        path = os.path.expanduser(path)

        if not os.path.exists(path):
            return output_error(msg('fail_path_doesnt_exist', path), cmd_pointer)
        cmd_pointer.settings['paths'][workspace_name] = path
    spinner.start('Creating workspace')
    sleep(0.5)  # Ensure the spinner is displayed for at least a moment.
    cmd_pointer.settings['workspaces'].append(workspace_name)
    cmd_pointer.settings['workspace'] = workspace_name
    cmd_pointer.histfile = os.path.expanduser(
        cmd_pointer.workspace_path(cmd_pointer.settings['workspace'].upper()) + '/.cmd_history'
    )

    error_creating_dir = False
    error_other = False
    try:
        workspace_name = cmd_pointer.settings['workspace'].upper()
        dir_path = os.path.expanduser(cmd_pointer.workspace_path(workspace_name))
        if not os.path.exists(dir_path):
            os.mkdir(dir_path)
        else:
            # This currently happens when you remove a workspace and then try to recreate it.
            # @Phil - we probably should move or archive the workspace folder when removing the workspace.
            os.chdir(dir_path)
            error_creating_dir = msg('warn_workspace_folder_already_exists', workspace_name)
        # Main and session registry writes
        write_registry(cmd_pointer.settings, cmd_pointer, True)
        write_registry(cmd_pointer.settings, cmd_pointer)
        
        readline.clear_history()
        readline.write_history_file(cmd_pointer.histfile)
        # raise ValueError('This is a test error.\n') @Phil this causes the app to break permamenently.
    except BaseException as err:
        error_other = msg('err_workspace_create', err, split=True)

    # Show success/errror message.
    if cmd_pointer.api_mode is False:
        if error_other:
            spinner.fail(output_error(error_other, cmd_pointer, return_val=True, jup_return_format='plain'))
            spinner.start()
            spinner.stop()
            return
        else:
            add_line = bool(error_creating_dir)
            error_creating_dir = error_creating_dir + '\n' if error_creating_dir else ''
            spinner.succeed(
                output_text(
                    msg('success_workspace_create', workspace_name, error_creating_dir),
                    cmd_pointer,
                    return_val=True,
                    jup_return_format='plain'
                )
            )
            spinner.start()
            spinner.stop()
            if add_line:
                print('')
            return
