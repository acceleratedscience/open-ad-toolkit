import os
from openad.helpers.files import open_file


# Todo: make list_files use get_workspace_files, has duplicate functionality now.
def fs_get_workspace_files(cmd_pointer, path=""):
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


def fs_get_file(cmd_pointer, path):
    """
    Read a file or directory from the workspace.
    """

    # Compile path
    workspace_path = cmd_pointer.workspace_path(cmd_pointer.settings["workspace"])
    file_path = workspace_path + "/" + path

    # Read file
    data, err_code = open_file(file_path, return_err="code", as_string=True)

    # if data:
    #     print(data)
    #     print(len(data))
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
