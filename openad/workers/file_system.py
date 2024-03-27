import os
import time
import json
import shutil
import asyncio
from openad.helpers.files import open_file, file_stats
from openad.helpers.timeit import timeit
from openad.gui.api.molecules_api import (
    MoleculesApi,
    get_molset_mols,
    create_molset_response,
    smiles_file2molset,
    sdf2molset,
    index_molset_file_async,
)


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
            file_path = path + "/" + filename if path else filename
            items[index] = fs_compile_filedir_obj(cmd_pointer, file_path)

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


def fs_compile_filedir_obj(cmd_pointer, path):
    """
    Compile universal file object that cen be parsed by the frontend.

    Used by the file browser (FileBroswer.vue) and the file viewers (FileStore).
    The file content is added later, under the "data" key - see fs_get_file_data().

    Parameters:
    -----------
    cmd_pointer: object
        The command pointer object, used to fetch the workspace path.
    path: str
        The path of the file relative to the workspace, including the filename.

    """
    workspace_path = cmd_pointer.workspace_path(cmd_pointer.settings["workspace"])
    path_absolute = os.path.join(workspace_path, path)
    filename = os.path.basename(path)

    # Get file exists or error code.
    f_stats, err_code = file_stats(path_absolute)

    # No file/dir found
    if err_code:
        return {
            "_meta": {"errCode": err_code},
            "filename": filename,
            "path": path,
            "pathAbsolute": path_absolute,
        }

    # File
    if os.path.isfile(path_absolute):
        size = f_stats.st_size
        time_edited = f_stats.st_mtime * 1000  # Convert to milliseconds for JS.
        time_created = f_stats.st_ctime * 1000  # Convert to milliseconds for JS.
        ext = _get_file_ext(filename)
        ext2 = _get_file_ext2(filename)
        file_type = _get_file_type(ext, ext2)
        return {
            "_meta": {
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

    # Directory
    elif os.path.isdir(path_absolute):
        return {
            "_meta": {
                "fileType": "dir",
            },
            "filename": filename,
            "path": path,
            "pathAbsolute": path_absolute,
        }


def fs_attach_file_data(cmd_pointer, file_obj, query=None):
    """
    Read the content of a file and attach it to the file object.

    This content will then be attached to the "data" key of the file object
    to be consumed by the frontend. For most file types, this is just the
    raw text content of the file, but for certain files like a molset, this
    will be a parsed object that includes additional data like pagination etc.

    This entry point is only used for opening files, once a molset (or potentially
    other editable file formats later) is opened, further querying and editing
    is handled by its own API endpoint - i.e. query_molset()

    Parameters:
    -----------
    path: str
        The path of the file relative to the workspace, including the filename.
    query: dict
        The query object, used to filter and sort the molset, and possible other
        file formats in the future.
    """

    path_absolute = file_obj["pathAbsolute"]
    file_type = file_obj["_meta"]["fileType"]
    ext = file_obj["_meta"]["ext"]

    # Molset --> Load molset object with first page data
    if file_type == "molset":

        # Step 1: Load or assemble the molset.
        # - - -

        # From molset JSON file, no formatting required.
        if ext == "json":
            molset, err_code = get_molset_mols(path_absolute)

        # From SMILES file
        elif ext == "smi":
            molset, err_code = smiles_file2molset(path_absolute)

        # From SDF file
        elif ext == "sdf":
            molset, err_code = sdf2molset(path_absolute)

        # Step 2: Store a working copy of the molset in the cache.
        # - - -

        if molset:
            cache_id = str(int(time.time() * 1000))
            cache_path = fs_assemble_cache_path(cmd_pointer, "molset", cache_id)

            # For JSON files, we can simply copy the original file (fast).
            if ext == "json":
                # timeit("copy_file")
                shutil.copy(path_absolute, cache_path)
                # timeit("copy_file", True)

                # Add indices to molecules in our working copy,
                # without blocking the thread.
                # timeit("index_wc")
                index_molset_file_async(cache_path)
                # timeit("index_wc", True)

            # For other formats, we need to write the molset object to disk (slow).
            else:
                # timeit("write_cache")
                with open(cache_path, "w", encoding="utf-8") as f:
                    json.dump(molset, f)
                # timeit("write_cache", True)

            # Step 3: Create the response object.
            # - - -

            # Filter, sort & paginate the molset, wrap it into
            # a response object and add the cache_id so further
            # operations can be performed on the working copy.
            data = create_molset_response(molset, query, cache_id)

        else:
            data = None

    # Everything else --> Load file content
    else:
        # Read file's content
        data, err_code = open_file(path_absolute, return_err="code")

    # Attach file content or error code
    if err_code:
        file_obj["_meta"]["errCode"] = err_code
    else:
        file_obj["data"] = data

    return file_obj


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
    # Single molecule files
    if ext in ["mol", "molecule", "pdb", "cif", "xyz", "mol2", "mmcif", "cml", "inchi"]:
        return "mol"

    # Molecule set files
    if ext in ["sdf", "smi"]:
        return "molset"

    # JSON files --> parse secondary extension
    elif ext in ["json", "cjson"]:
        # Molecule
        if ext2 == "mol":
            return "mol"
        # Molecule set
        elif ext2 == "molset":
            return "molset"
        # JSON files
        else:
            return "json"

    # Data files
    elif ext in ["csv"]:
        return "data"

    # Text files
    elif ext in ["txt", "md", "yaml", "yml"]:
        return "text"

    # HTML files
    elif ext in ["html", "htm"]:
        return "html"

    # Image formats
    elif ext in ["jpg", "jpeg", "png", "gif", "bmp", "webp"]:
        return "img"

    # Video formats
    elif ext in ["mp4", "avi", "mov", "mkv", "webm"]:
        return "vid"

    # Individually recognized file formats (have their own icon)
    elif ext in ["xml", "pdf", "svg", "run", "rxn", "md"]:
        return ext

    # # Yaml files
    # elif ext in ["yaml", "yml"]:
    #     return "yaml"

    else:
        # Unrecognized file formats
        return "unk"


# Compile the file path to a cached working copy of a file.
def fs_assemble_cache_path(cmd_pointer, file_type, cache_id):
    """
    Compile the file path to a cached working copy of a file.

    Parameters:
    -----------
    cmd_pointer: object
        The command pointer object, used to fetch the workspace path.
    file_type: 'molset'
        The type of file, used to name the cache file. For now only molset.
    """
    workspace_path = cmd_pointer.workspace_path()
    return f"{workspace_path}/.cache/{file_type}-{cache_id}.json"
