"""Handles file System interactions"""
#!/usr/local/opt/python@3.9/bin/python3.9
#

import os
import shutil
from datetime import datetime

# Expand user path: ~/ --> ../
from pathlib import PosixPath

# Helpers
from openad.helpers.general import confirm_prompt
from openad.helpers.output import msg, output_text, output_error, output_success, output_table

# Globals

from openad.app.global_var_lib import _date_format


# Importing our own plugins.
# This is temporary until every plugin is available as a public pypi package.


def list_files(cmd_pointer, parser):
    """List all Files in a Workspace"""
    workspace_path = cmd_pointer.workspace_path(cmd_pointer.settings["workspace"])

    files = []
    table_headers = ("File Name", "Size", "Last Edited")

    # Directory from which we're reading files.
    path = workspace_path

    # Get list of all files only in the given directory.
    def fun(x):
        return os.path.isfile(os.path.join(path, x))

    files_list = filter(fun, os.listdir(path))

    # Create a list of tuples with file info: (filename, size, time)
    files_data = [
        (os.path.basename(f), os.stat(os.path.join(path, f)).st_size, os.stat(os.path.join(path, f)).st_atime)
        for f in files_list
    ]

    # Check if there are any non-hidden files in the workspace.
    non_hidden_files = False
    for file in files_data:
        if not file[0].startswith("."):
            non_hidden_files = True
            break

    # Display message when no files are found.
    # len(files_data) == 1
    if not non_hidden_files:
        workspace_name = cmd_pointer.settings["workspace"].upper()
        return output_text(msg("no_workspace_files", workspace_name), cmd_pointer, pad=1)

    # Assemble table data.
    for name, size, timestamp in sorted(files_data, key=lambda x: x[1], reverse=True):
        if name.startswith("."):
            # For now we're jumping over hidden files, though
            # I would like to add an option to display them.
            # Probably `list all files` - moenen
            continue

        if size < (1024 * 1024) / 10:
            size = f"{round(size / 1024, 2)} kB"
        else:
            size = f"{round(size / (1024 * 1024), 2)} MB"
        timestamp = datetime.fromtimestamp(timestamp)
        timestamp = timestamp.strftime(_date_format)
        result = [name, size, timestamp]
        files.append(result)

    # Display/return table.
    return output_table(files, cmd_pointer, headers=table_headers)


# External path to workspace path
def import_file(cmd_pointer, parser):
    """Import a file from thefiles system external to Workspaces"""
    # Reset working directory as it can have changed.
    # os.chdir(_repo_dir)

    workspace_path = cmd_pointer.workspace_path(cmd_pointer.settings["workspace"])
    workspace_name = cmd_pointer.settings["workspace"].upper()
    source_file = parser["source"]
    dest_file = parser["destination"]

    # Expand user path: ~/ --> ../
    # from pathlib import PosixPath # Trash
    # source_file = PosixPath(source_file).expanduser().resolve() # Trash
    source_file = os.path.expanduser(source_file)

    if not os.path.exists(source_file):
        # Source does not exist
        return output_error(msg("fail_file_doesnt_exist", source_file), cmd_pointer)
    elif os.path.exists(workspace_path + "/" + dest_file):
        # Destination already exists
        if not confirm_prompt("Destination file already exists. Overwrite?"):
            return output_error(msg("abort"), cmd_pointer)
    try:
        # Success
        # shutil.copyfile(PosixPath(source_file).expanduser().resolve(), path + '/' + dest_file) # Trash
        shutil.copyfile(os.path.expanduser(source_file), workspace_path + "/" + dest_file)
        return output_success(msg("success_import", source_file, workspace_name), cmd_pointer)
    except Exception as err:
        # Failure
        return output_error(msg("err_import", err, split=True), cmd_pointer)


# Workspace path to external path
def export_file(cmd_pointer, parser):
    """Exports a workspace file to the rechable filesystem"""
    # Reset working directory as it can have changed.
    # os.chdir(_repo_dir)

    workspace = cmd_pointer.workspace_path(cmd_pointer.settings["workspace"])
    source_file = parser["source"]
    dest_file = parser["destination"]
    workspace_name = cmd_pointer.settings["workspace"].upper()

    dest_file = PosixPath(dest_file).expanduser().resolve()

    if not os.path.exists(workspace + "/" + source_file):
        # Source does not exist
        return output_error(msg("fail_file_doesnt_exist", workspace + "/" + source_file), cmd_pointer)

    elif os.path.exists(dest_file) is True:
        # Destination already exists
        if not confirm_prompt("Destination file already exists. Overwrite?"):
            return output_error(msg("abort"), cmd_pointer)
    try:
        # Success
        shutil.copyfile(workspace + "/" + source_file, dest_file)
        return output_success(msg("success_export", source_file, workspace_name, dest_file), cmd_pointer)
    except Exception as err:
        # Failure
        return output_error(msg("err_export", err, split=True), cmd_pointer)


# Workspace path to workspace name
def copy_file(cmd_pointer, parser):
    """copy a file betqeen workspaces"""
    # Reset working directory as it can have changed.
    # os.chdir(_repo_dir)

    source_file = parser["source"]
    source_file_path = cmd_pointer.workspace_path(cmd_pointer.settings["workspace"]) + "/" + source_file
    dest_file_path = cmd_pointer.workspace_path(parser["destination"]) + "/" + source_file
    source_workspace_name = cmd_pointer.settings["workspace"].upper()
    dest_workspace_name = parser["destination"].upper()

    if not os.path.exists(source_file_path):
        # Source does not exist
        return output_error(msg("fail_file_doesnt_exist", source_file_path), cmd_pointer)
    elif (
        parser["destination"].upper() != source_workspace_name
        and dest_workspace_name not in cmd_pointer.settings["workspaces"]
    ):
        # Invalid destination
        return output_error(msg("invalid_workpace_destination", parser["destination"].upper()), cmd_pointer)
    elif os.path.exists(dest_file_path) is True:
        # Destination already exists
        if not confirm_prompt("Destination file already exists. Overwrite?"):
            return output_error(msg("abort"), cmd_pointer)
    try:
        # Success
        shutil.copyfile(source_file_path, dest_file_path)
        return output_success(msg("success_copy", source_file, source_workspace_name, dest_workspace_name), cmd_pointer)
    except Exception as err:
        # Failure
        return output_error(msg("err_copy", err, split=True), cmd_pointer)


# Workspace path
def remove_file(cmd_pointer, parser):
    """remove a file from a workspace"""
    workspace = cmd_pointer.workspace_path(cmd_pointer.settings["workspace"])
    file_name = parser["file"]
    file_path = workspace + "/" + file_name
    workspace_name = cmd_pointer.settings["workspace"].upper()

    if not os.path.exists(file_path):
        # Source does not exist
        return output_error(msg("fail_file_doesnt_exist", file_path), cmd_pointer)
    if not confirm_prompt("Are you sure? This cannot be undone."):
        # Confirm prompt
        return output_error(msg("abort"), cmd_pointer)
    try:
        # Success
        os.remove(file_path)
        return output_success(msg("success_delete", file_name, workspace_name), cmd_pointer)
    except Exception as err:
        # Failure
        return output_error(msg("err_delete", err, split=True), cmd_pointer)
