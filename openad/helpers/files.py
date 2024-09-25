import os
import json
import ijson
import pandas as pd
from io import StringIO
from openad.helpers.output import output_error
from openad.helpers.output_msgs import msg
from openad.helpers.json_decimal_encoder import DecimalEncoder


# Standardized file opener.
# Detects file type and returns appropriate data format:
# - JSON/CJSON: dict or list of dicts
# - CSV: pandas dataframe
# - Other: string
def open_file(file_path, return_err=False, dumb=False, page=None, _stats_only=False):
    """
    Takes care of all boilerplate file opening code.

    Parameters:
    -----------
    file_path: str
        The path of the file to open.
    return_err: bool | 'code'
        Return a tuple (data, err_msg) instead of printing the error.
    dumb: bool
        To be deleted after testing, we don't want to return dataframes ever
    _stats_only: bool
        Only check if the file can be opened, don't return the data.
    """
    ext = file_path.split(".")[-1].lower()
    data = None
    err_msg = None
    err_code = None
    try:
        # Check if the file can be opened. This should only be invoked via file_stats.
        if _stats_only:
            data = os.stat(file_path)

        # Return file content.
        else:
            # Large files - load in chunks
            if _is_large_file(file_path):
                chunk_size = 10**3  # 1000 - The number of rows (csv) or objects (json) to load at a time.
                start = chunk_size * (page - 1) if page else 0

                # Parse JSON
                if ext == "json" or ext == "cjson":
                    data = get_chunk_json(file_path, start, chunk_size, dumb)

                # Parse CSV
                elif ext == "csv":
                    data = get_chunk_text(file_path, start, chunk_size, dumb, is_csv=True)

                # Parse any text file
                else:
                    data = get_chunk_text(file_path, start, chunk_size, dumb)

            # Regular files - load entirely
            else:
                with open(file_path, "r", encoding="utf-8") as f:
                    # Return string
                    if dumb:
                        data = f.read()

                    # Parse JSON as dict
                    elif ext == "json" or ext == "cjson":
                        data = json.load(f)

                    # # Parse CSV as dataframe
                    # elif ext == "csv":
                    #     data = pd.read_csv(f)

                    # Return string
                    else:
                        data = f.read()

        # Return data
        if return_err:
            return data, None
        else:
            return data
    except FileNotFoundError:
        err_msg = msg("err_file_not_found", file_path)
        err_code = "not_found"
    except PermissionError:
        err_msg = msg("err_file_no_permission_read", file_path)
        err_code = "no_permission"
    except IsADirectoryError:
        err_msg = msg("err_file_is_dir", file_path)
        err_code = "is_dir"
    except UnicodeDecodeError:
        err_msg = msg("err_decode", file_path)
        err_code = "decode"
    except IOError as err:
        err_msg = msg("err_io", file_path, err.strerror)
        err_code = "io"
    except BaseException as err:
        err_msg = msg("err_unknown", err)
        err_code = "unknown"
    # Note: if ever any new error codes are added here,
    # they should be mirrored in the openad-gui's FileStore.

    # Return error
    if return_err:
        if return_err == "code":
            return None, err_code
        else:
            return None, err_msg

    # Display error
    else:
        output_error(err_msg)
        return None


# Disabled until we have time to handle this properly.
# The idea is that when you try to open a mega file, we
# let you open it in chunks of 10k (or something) items
# at a time, but right now we just truncate the items
# without any warning which is bad.
def _is_large_file(file_path):
    return False
    file_size = os.path.getsize(file_path)
    return bool(file_size > 2 * 1024 * 1024)


# Get a file's stats (os.stat) or error code if invalid.
def file_stats(file_path, return_err="code"):
    return open_file(file_path, return_err, _stats_only=True)


# Read only a chunk of a large JSON file.
# This assumes the JSON consists of a list of objects.
def get_chunk_json(file_path, start, chunk_size, as_string=False):
    with open(file_path, "r", encoding="utf-8") as f:
        objects = ijson.items(f, "item")
        chunk = []
        for i, obj in enumerate(objects):
            if i < start:
                continue
            if i >= start + chunk_size:
                break
            chunk.append(obj)

        # If the target file is a single object, we can't iterate.
        # Instead, ijson.items will return an empty list, in which
        # case we read the entire file and return it as a single chunk.
        if chunk == []:
            f.seek(0)
            chunk = [json.loads(f.read())]

    if as_string:
        return json.dumps(chunk, cls=DecimalEncoder)
    else:
        return chunk


# Read only a chunk of a large text file.
def get_chunk_text(file_path, start, chunk_size, as_string=False, is_csv=False):
    with open(file_path, "r", encoding="utf-8") as f:
        chunk = []

        for i, line in enumerate(f):
            # Always include the first line as header.
            if i == 0:
                chunk.append(line.strip())
                continue
            # Then +1 to not count the header line.
            if i < start + 1:
                continue
            if i >= start + 1 + chunk_size:
                break
            chunk.append(line.strip())

    # Return CSV as pandas dataframe.
    if is_csv and not as_string:
        return pd.read_csv(StringIO("\n".join(chunk)))
    else:
        return "\n".join(chunk)


# TO DO: implement this
def count_items_in_json(file_path):
    with open(file_path, "r", encoding="utf-8") as f:
        objects = ijson.items(f, "item")
        return sum(1 for _ in objects)


# # Read an entire JSON file.
# def get_data_json(file_path):
#     with open(file_path, "r", encoding="utf-8") as f:
#         data = json.load(f)
#     return data


# # Read an entire CSV file. -- as dict
# def get_data_csv(file_path):
#     with open(file_path, "r", encoding="utf-8") as f:
#         data = csv.DictReader(f)
#         data = list(data)
#     return data


# # Read only a chunk of a large CSV file. -- as dicttionaries
# def get_chunk_csv(file_path, start, chunk_size):
#     with open(file_path, "r", encoding="utf-8") as f:
#         reader = csv.DictReader(f)
#         chunk = []
#         for i, row in enumerate(reader):
#             if i < start:
#                 continue
#             if i >= start + chunk_size:
#                 break
#             chunk.append(row)
#         return chunk


# Standardized file writer.
def write_file(file_path, data, return_err=False):
    err_msg = None
    try:
        with open(file_path, "w", encoding="utf-8") as f:
            f.write(data)

            # Return success
            if return_err:
                return True, None
            else:
                return True
    except FileNotFoundError:
        err_msg = msg("err_file_not_found", file_path)
    except PermissionError:
        err_msg = msg("err_file_no_permission_write", file_path)
    except IsADirectoryError:
        err_msg = msg("err_file_is_dir", file_path)
    except UnicodeDecodeError:
        err_msg = msg("err_decode", file_path)
    except IOError as err:
        err_msg = msg("err_io", file_path, err.strerror)
    except Exception as err:
        err_msg = msg("err_unknown", err)

    # Return error
    if return_err:
        return None, err_msg

    # Display error
    else:
        output_error(err_msg)
        return None


# Remove trash folder and its contents.
# Called on quit.
def empty_trash(cmd_pointer):
    trash_dir = f"{cmd_pointer.workspace_path()}/.trash"
    os.system(f"rm -rf '{trash_dir}'")
