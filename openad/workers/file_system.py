import os
from openad.helpers.files import open_file, file_stats
from openad.smols.smol_functions import create_molset_cache_file, get_molset_mols
from openad.gui.api.molecules_api import create_molset_response
from openad.smols.smol_transformers import smiles_path2molset, sdf_path2molset, mdl_path2smol
from openad.mmols.mmol_transformers import cif2mmol, pdb2mmol


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
        # "filesSystem": [],  # Filenames starting with __ #fileSystem A place to hide our own system files out of sight? Will probably implement later, do not delete
        "dirs": [],
        "dirsHidden": [],  # Dir names starting with .
        # "dirsSystem": [],  # Dir names starting with __ # See #fileSystem
    }

    # Organize file & directory names into dictionary.
    for filename in os.listdir(dir_path):
        is_hidden = filename.startswith(".")
        # is_system = filename.startswith("__")
        is_file = os.path.isfile(os.path.join(dir_path, filename))

        if is_file:
            if filename == ".DS_Store":
                continue
            elif is_hidden:
                level["filesHidden"].append(filename)
            # elif is_system:
            #     level["filesSystem"].append(filename) # See #fileSystem
            else:
                level["files"].append(filename)
        else:
            is_dir = os.path.isdir(os.path.join(dir_path, filename))
            if is_dir:
                if filename == "._openad":
                    continue
                if is_hidden:
                    level["dirsHidden"].append(filename)
                # elif is_system:
                #     level["dirsSystem"].append(filename) # See #fileSystem
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
    is handled by its own API endpoint - i.e. get_molset()

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

    # Molset files --> Load molset object with first page data
    if file_type in ["molset", "sdf", "smi"]:
        # Step 1: Load or assemble the molset.
        # - - -

        # From molset JSON file, no formatting required.
        if ext == "json":
            molset, err_code = get_molset_mols(path_absolute)

        # From SMILES file
        elif ext == "smi":
            molset, err_code = smiles_path2molset(path_absolute)

        # From SDF file
        elif ext == "sdf":
            molset, err_code = sdf_path2molset(path_absolute)

        if molset:
            # Step 2: Store a working copy of the molset in the cache.
            # - - -

            # For JSON files, we can simply copy the original file (fast).
            if ext == "json":
                cache_id = create_molset_cache_file(cmd_pointer, path_absolute=path_absolute)

            # All other cases, write file from memory.
            else:
                cache_id = create_molset_cache_file(cmd_pointer, molset=molset)

            # Step 3: Create the response object.
            # - - -

            # Filter, sort & paginate the molset, wrap it into
            # a response object and add the cache_id so further
            # operations can be performed on the working copy.
            data = create_molset_response(molset, query, cache_id)

        else:
            data = None

    # Molecule files --> convert to molecule JSON
    elif file_type in ["mdl", "pdb", "cif"]:
        # From MOL file
        if ext == "mol":
            data, err_code = mdl_path2smol(path_absolute)

        # From PDB file
        if ext == "pdb":
            data = pdb2mmol(pdb_path=path_absolute)
            err_code = None

        # From CIF file
        if ext == "cif":
            data = cif2mmol(cif_path=path_absolute)
            err_code = None

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
    """
    Get the file extension from a filename.
    """
    if filename.find(".") == -1:
        return ""
    else:
        return filename.split(".")[-1]


def _get_file_ext2(filename):
    """
    Get the secondary file extension from a filename.

    Secondary file extensions are used to indicate subformats, eg:
    - foobar.json --> JSON file
    - foobar.mol.json --> molecule JSON file
    - foobar.molset.json --> molecule set JSON file
    """
    has_ext2 = len(filename.split(".")) >= 3
    if has_ext2:
        parts = filename.split(".")
        parts.pop()
        return parts.pop() or None
    else:
        return None


def _get_file_type(ext, ext2):
    """
    Deduct the fileType from the file's primary and secondary extensions.

    The file type is used in the frontend to determine:
    - What icon to display
    - What viewer to use when opening the file

    In the frontend the fileType is mapped to its respective
    display name and module name via _map_FileType.

    Any changes here should also be reflected in the FileType TypeScript type.
    """
    # Small molecule files
    if ext in ["mol"]:  # Future support: "molecule", "pdb", "cif", "xyz", "mol2", "mmcif", "cml", "inchi"
        return "mdl"

    # Macromolecule files
    if ext in ["pdb"]:
        return "pdb"
    if ext in ["cif"]:
        return "cif"

    # Molecule set files
    if ext in ["smi"]:
        return "smi"

    # JSON files --> parse secondary extension
    elif ext in ["json", "cjson"]:
        # Small molecule
        if ext2 == "smol":
            return "smol"
        elif ext2 == "mol":  # backward compatibility for mol.json files
            return "smol"
        elif ext2 == "mmol":
            return "mmol"
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
    elif ext in ["sdf", "xml", "pdf", "svg", "run", "rxn", "md"]:
        return ext

    # # Yaml files
    # elif ext in ["yaml", "yml"]:
    #     return "yaml"

    else:
        # Unrecognized file formats
        return "unk"
