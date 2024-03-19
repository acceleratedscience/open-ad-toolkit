import os
from openad.helpers.files import open_file
from openad.gui.api.molecules_api import MoleculesApi


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
            "empty": False,
            "empty_hidden": False,
        },
        "dirname": "",
        "files": [],
        "filesHidden": [],  # Filenames starting with .
        # "filesSystem": [],  # Filenames starting with __ # Probably we can just use hidden for this.
        "dirs": [],
        "dirsHidden": [],  # Dir names starting with .
        # "dirsSystem": [],  # Dir names starting with __ # Probably we can just use hidden for this.
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
                level["filesHidden"].append(filename)
            elif is_system:
                level["filesSystem"].append(filename)
            else:
                level["files"].append(filename)
        else:
            is_dir = os.path.isdir(os.path.join(dir_path, filename))
            if is_dir:
                if is_hidden:
                    level["dirsHidden"].append(filename)
                elif is_system:
                    level["dirsSystem"].append(filename)
                else:
                    level["dirs"].append(filename)

    # Sort the lists
    level["files"].sort()
    level["filesHidden"].sort()
    level["dirs"].sort()
    level["dirsHidden"].sort()

    # Expand every dir & filename into a dictionary: {_meta, filename, path}
    for category, items in level.items():
        if category == "_meta":
            continue
        for index, filename in enumerate(items):
            is_dir = category in ["dirs", "dirsHidden"]
            file_path = path + "/" + filename if path else filename
            items[index] = _compile_file_data(cmd_pointer, file_path, filename, is_dir)

    #
    #

    # Attach workspace name
    workspace_name = cmd_pointer.settings["workspace"].upper()
    dir_name = path.split("/")[-1]
    level["dirname"] = workspace_name if not path else dir_name

    # Mark empty directories.
    if not level["files"] and not level["dirs"]:
        level["_meta"]["empty"] = True
    if level["_meta"]["empty"] and not level["filesHidden"] and not level["dirsHidden"]:
        level["_meta"]["empty_hidden"] = True

    return level


def fs_get_file(cmd_pointer, path):
    """
    Read a file or directory from the workspace.
    """

    # Filepath
    workspace_path = cmd_pointer.workspace_path(cmd_pointer.settings["workspace"])
    path_absolute = workspace_path + "/" + path
    filename = path.split("/")[-1]

    # Check if path is a file or a directory
    file_or_dir = _check_path(path_absolute)

    # Not found
    if file_or_dir is None:
        return {"_meta": {}}

    # Directory
    elif file_or_dir == "dir":
        return {
            "_meta": {
                "fileType": "dir",  # "dir_hidden" if filename[0] == "." else "dir",
                "path": path,
            }
        }

    # File
    elif file_or_dir == "file":
        file = _compile_file_data(cmd_pointer, path, filename)

        # Molset --> Load molset object with first page data
        if file["_meta"]["fileType"] == "molset":
            molecules_api = MoleculesApi(cmd_pointer)
            data = molecules_api.get_molset()
            file["data"] = data

        # Everything else --> Load file content
        else:
            # Read file's content
            data, err_code = open_file(path_absolute, return_err="code", dumb=False)  # %%

            # Attach file content or error code
            if err_code:
                file["_meta"]["errCode"] = err_code
            else:
                file["data"] = data

        return file


def _check_path(path):
    """
    Check if a path is a file or a directory.
    """
    if os.path.isdir(path):
        return "dir"
    elif os.path.isfile(path):
        return "file"
    else:
        return None


def _compile_file_data(cmd_pointer, path, filename, is_dir=False):
    """
    Compile universal file object that will be parsed by the frontend.
    Used by the file browser (FileBroswer.vue) and the file viewers (FileStore).

    Parameters:
    -----------
    cmd_pointer: object
        The command pointer object, used to fetch the workspace path.
    path: str
        The path of the file relative to the workspace, including the filename.
    filename: str
        The name of the file.

    """
    workspace_path = cmd_pointer.workspace_path(cmd_pointer.settings["workspace"])
    path_absolute = workspace_path + "/" + path
    size = os.stat(path_absolute).st_size
    time_edited = os.stat(path_absolute).st_mtime * 1000  # Convert to milliseconds for JS.
    time_created = os.stat(path_absolute).st_ctime * 1000  # Convert to milliseconds for JS.
    ext = "" if is_dir else _get_file_ext(filename)
    ext2 = "" if is_dir else _get_file_ext2(filename)
    file_type = "dir" if is_dir else _get_file_type(ext, ext2)

    # Dict structure for one file/dir.
    return {
        "_meta": {
            # "name": filename,
            "fileType": file_type,
            "ext": ext,
            "ext2": ext2,  # Secondary file extension, eg. foobar.mol.json --> mol
            "size": size,
            "timeCreated": time_created,
            "timeEdited": time_edited,
            "errCode": None,
        },
        "data": None,  # Just for reference, this is added when opening file.
        "filename": filename,
        "path": path,  # Relative to workspace.
        "pathAbsolute": path_absolute,  # Absolute path
    }


def _get_file_ext(filename):
    if filename.find(".") == -1:
        return ""
    else:
        return filename.split(".")[-1]


# Secondary file extension, eg. foobar.mol.json --> mol
def _get_file_ext2(filename):
    has_ext2 = len(filename.split(".")) >= 3
    if has_ext2:
        parts = filename.split(".")
        parts.pop()
        return parts.pop() or None
    else:
        return None


def _get_file_type(ext, ext2):
    if ext in ["sdf", "mol", "molecule", "pdb", "cif", "xyz", "mol2", "mmcif", "cml", "smiles", "inchi"]:
        # Molecule formats
        return "mol"
    elif ext in ["csv"]:
        # Data formats
        return "data"
    elif ext in ["json", "cjson"]:
        if ext2 == "mol":
            # Molecule
            return "mol"
        elif ext2 == "molset":
            # Molecule set
            return "molset"
        else:
            # JSON files
            return "json"
    elif ext in ["txt", "md", "yaml", "yml"]:
        # Text formats
        return "text"
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
        return "unk"
