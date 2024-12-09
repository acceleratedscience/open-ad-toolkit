"""Handles file System interactions"""

#!/usr/local/opt/python@3.9/bin/python3.9
#

import os
import shutil
from datetime import datetime

# Expand user path: ~/ --> ../
from pathlib import PosixPath

# Workers
from openad.workers.file_system import fs_get_workspace_files

# Helpers
from openad.helpers.general import confirm_prompt
from openad.helpers.output import output_text, output_error, output_success, output_table, strip_tags
from openad.helpers.output_msgs import msg


# Globals
from openad.app.global_var_lib import _date_format


# Importing our own plugins.
# This is temporary until every plugin is available as a public pypi package.


def list_files(cmd_pointer, parser):
    import pprint

    path = parser["path"] if "path" in parser else ""
    data = fs_get_workspace_files(cmd_pointer, path)
    space = [""] if data["dirs"] else []
    files = data["dirs"] + space + data["files"]
    # pprint.pprint(files)
    table = []
    table_headers = ("File Name", "Size", "Last Edited")
    for file in files:
        # Insert space
        if not file:
            table.append(("-", "-", "-"))
            continue

        filetype = file["_meta"]["fileType"]
        filename = file["filename"] + ("/" if filetype == "dir" else "")
        size = file["_meta"]["size"] if "size" in file["_meta"] else None
        timestamp = file["_meta"]["timeEdited"] if "timeEdited" in file["_meta"] else None

        if filename.startswith("."):
            # For now we're jumping over hidden files, though
            # I would like to add an option to display them.
            # Probably `list all files` - moenen
            continue

        if size:
            if size < (1024 * 1024) / 10:
                size = f"{round(size / 1024, 2)} kB"
            else:
                size = f"{round(size / (1024 * 1024), 2)} MB"

        if timestamp:
            timestamp = datetime.fromtimestamp(timestamp / 1000)
            timestamp = timestamp.strftime(_date_format)

        result = (filename, size, timestamp)
        table.append(result)

    # return "OK"
    return output_table(table, is_data=False, headers=table_headers, colalign=("left", "right", "left"))


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
        return output_error(msg("err_file_doesnt_exist", source_file))
    elif os.path.exists(workspace_path + "/" + dest_file):
        # Destination already exists
        if not confirm_prompt("Destination file already exists. Overwrite?"):
            return output_error(msg("abort"))

    try:
        # Success
        # shutil.copyfile(PosixPath(source_file).expanduser().resolve(), path + '/' + dest_file) # Trash
        if os.path.isfile(source_file):
            shutil.copyfile(source_file, workspace_path + "/" + dest_file)
        else:
            # @later: Change language to reflect dir instead of file
            shutil.copytree(source_file, workspace_path + "/" + dest_file)
        return output_success(msg("success_import", source_file, workspace_name))
    except Exception as err:
        # Failure
        return output_error(msg("err_import", err))


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
        return output_error(msg("err_file_doesnt_exist", workspace + "/" + source_file))

    elif os.path.exists(dest_file) is True:
        # Destination already exists
        if not confirm_prompt("Destination file already exists. Overwrite?"):
            return output_error(msg("abort"))
    try:
        # Success
        shutil.copyfile(workspace + "/" + source_file, dest_file)
        return output_success(msg("success_export", source_file, workspace_name, dest_file))
    except Exception as err:
        # Failure
        return output_error(msg("err_export", err))


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
        return output_error(msg("err_file_doesnt_exist", source_file_path))
    elif (
        parser["destination"].upper() != source_workspace_name
        and dest_workspace_name not in cmd_pointer.settings["workspaces"]
    ):
        # Invalid destination
        return output_error(msg("invalid_workpace_destination", parser["destination"].upper()))
    elif os.path.exists(dest_file_path) is True:
        # Destination already exists
        if not confirm_prompt("Destination file already exists. Overwrite?"):
            return output_error(msg("abort"))
    try:
        # Success
        shutil.copyfile(source_file_path, dest_file_path)
        return output_success(msg("success_copy", source_file, source_workspace_name, dest_workspace_name))
    except Exception as err:
        # Failure
        return output_error(msg("err_copy", err))


# Workspace path
def remove_file(cmd_pointer, parser):
    """remove a file from a workspace"""
    workspace = cmd_pointer.workspace_path(cmd_pointer.settings["workspace"])
    file_name = parser["file"]
    file_path = workspace + "/" + file_name
    workspace_name = cmd_pointer.settings["workspace"].upper()

    if not os.path.exists(file_path):
        # Source does not exist
        return output_error(msg("err_file_doesnt_exist", file_path))
    if not confirm_prompt("Are you sure? This cannot be undone."):
        # Confirm prompt
        return output_error(msg("abort"))
    try:
        # Success
        os.remove(file_path)
        return output_success(msg("success_delete", file_name, workspace_name))
    except Exception as err:
        # Failure
        return output_error(msg("err_delete", err))


def open_file(cmd_pointer, parser):
    from openad.gui.gui_launcher import gui_init

    path = "~/" + parser["file"]
    gui_init(cmd_pointer, path)
