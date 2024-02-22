"""Handles file System interactions"""

#!/usr/local/opt/python@3.9/bin/python3.9
#

import os
import shutil
from datetime import datetime

# Expand user path: ~/ --> ../
from pathlib import PosixPath

# Helpers
from openad.helpers.general import confirm_prompt, open_file
from openad.helpers.output import output_text, output_error, output_success, output_table, strip_tags
from openad.helpers.output_msgs import msg


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
        return output_text(msg("no_workspace_files", workspace_name), pad=1)

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
    return output_table(files, is_data=False, headers=table_headers)


# Todo: make list_files use get_workspace_files, has duplicate functionality now.
def get_workspace_files(cmd_pointer, path=""):
    """
    Return your active workspace's content as a JSON object.
    """

    # Get workspace path.
    workspace_path = cmd_pointer.workspace_path(cmd_pointer.settings["workspace"])
    dir_path = workspace_path + "/" + path

    # Dict structure for one level.
    level = {
        "_meta": {
            "name": "",
            "empty": False,
            "empty_hidden": False,
        },
        "files": [],
        "files_hidden": [],  # Filenames starting with .
        # "files_system": [],  # Filenames starting with __  # Probably we can just use hidden for this.
        "dirs": [],
        "dirs_hidden": [],  # Dir names starting with .
        # "dirs_system": [],  # Dir names starting with __ # Probably we can just use hidden for this.
    }

    # Organize file & directory names into dictionary.
    for filename in os.listdir(dir_path):
        is_hidden = filename.startswith(".")
        is_system = filename.startswith("__")
        is_file = os.path.isfile(os.path.join(dir_path, filename))

        if is_file:
            if filename == ".DS_Store":
                continue
            elif is_hidden:
                level["files_hidden"].append(filename)
            elif is_system:
                level["files_system"].append(filename)
            else:
                level["files"].append(filename)
        else:
            is_dir = os.path.isdir(os.path.join(dir_path, filename))
            if is_dir:
                if is_hidden:
                    level["dirs_hidden"].append(filename)
                elif is_system:
                    level["dirs_system"].append(filename)
                else:
                    level["dirs"].append(filename)

    # Sort the lists
    level["files"].sort()
    level["files_hidden"].sort()
    level["dirs"].sort()
    level["dirs_hidden"].sort()

    # Expand every dir & filename into a dictionary: {_meta, filename, path}
    for category, items in level.items():
        if category == "_meta":
            continue
        for index, filename in enumerate(items):
            path_full = os.path.join(dir_path, filename)
            path_relative = path + ("/" if path else "") + filename
            size = os.stat(path_full).st_size
            time_edited = os.stat(path_full).st_mtime * 1000
            time_created = os.stat(path_full).st_ctime * 1000
            file_ext = _get_file_ext(category, filename)
            file_type = _get_file_type(category, file_ext)

            # Dict structure for one file/dir.
            items[index] = {
                "_meta": {
                    "name": filename,
                    "size": size,
                    "time_edited": time_edited,
                    "time_created": time_created,
                    "type": file_type,
                    "ext": file_ext,
                },
                "filename": filename,
                "path": path_relative,
            }

    #
    #

    # Attach workspace name
    workspace_name = cmd_pointer.settings["workspace"].upper()
    dir_name = path.split("/")[-1]
    level["_meta"]["name"] = workspace_name if not path else dir_name

    # Mark empty directories.
    if not level["files"] and not level["dirs"]:
        level["_meta"]["empty"] = True
    if level["_meta"]["empty"] and not level["files_hidden"] and not level["dirs_hidden"]:
        level["_meta"]["empty_hidden"] = True

    return level


def _get_file_ext(category, filename):
    if category in ["dirs", "dirs_hidden"]:
        return ""
    elif filename.find(".") == -1:
        return ""
    else:
        return filename.split(".")[-1]


def _get_file_type(category, ext):
    if category in ["dirs", "dirs_hidden"]:
        # Directories
        return "dir"
    if ext in ["sdf", "mol", "molecule", "pdb", "cif", "xyz", "mol2", "mmcif", "cml", "smiles", "inchi"]:
        # Molecule formats
        return "mol"
    elif ext in ["csv"]:
        # Data formats
        return "data"
    elif ext in ["json", "cjson"]:
        # JSON files
        return "json"
    elif ext in ["txt", "md", "yaml", "yml"]:
        # Text formats
        return "txt"
    elif ext in ["xml", "pdf", "svg", "run", "rxn", "mod"]:
        # Individually recognized file formats (have their own icon)
        return ext
    elif ext in ["html", "htm"]:
        # HTML files
        return "html"
    elif ext in ["jpg", "jpeg", "png", "gif", "bmp", "webp"]:
        # Image formats
        return "img"
    elif ext in ["mp4", "avi", "mov", "mkv", "webm"]:
        # Video formats
        return "vid"
    # elif ext in ["yaml", "yml"]:
    #     # Yaml files
    #     return "yaml"
    else:
        # Unrecognized file formats
        return "doc"


def get_file(cmd_pointer, path):
    """
    Read a file or directory from the workspace.
    """

    # Compile path
    workspace_path = cmd_pointer.workspace_path(cmd_pointer.settings["workspace"])
    file_path = workspace_path + "/" + path

    # Read file
    data, err_code = open_file(file_path, return_err="code", raw=True)
    if err_code is None:
        # File content
        return {
            "data": data,
            "path": path,
            "pathFull": file_path,
            "isDir": False,
            "errCode": None,
        }
    elif err_code == "is_dir":
        # Directory
        return {
            "data": None,
            "path": path,
            "pathFull": file_path,
            "isDir": True,
            "errCode": None,
        }
    else:
        # File error
        return {
            "data": None,
            "path": path,
            "pathFull": file_path,
            "isDir": False,
            "errCode": err_code,
        }


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
