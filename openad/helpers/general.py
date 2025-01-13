import os
import re
import sys
import time
import shutil
import getpass
import readline
from datetime import datetime
from IPython.display import clear_output
from openad.helpers.output import output_text, output_error
from openad.helpers.output_msgs import msg
from openad.plugins.style_parser import style


# Refreshes the command prompt when in the shell.
def refresh_prompt(settings):
    if settings["context"] is not None:
        # prompt = ' \u001b[7m ' + settings['context'] + ' \u001b[0m '  # Reverse & reset
        prompt = settings["context"] + "->"  # Reverse & reset
    else:
        prompt = "OpenAD:"
    if settings["workspace"] is not None:
        prompt = prompt + settings["workspace"]
    prompt = prompt + " >>  "
    return prompt


# Allows us to detect Jupyter before GLOBAL_SETTINGS["display"] is set.
def is_notebook_mode():
    """Return True if we are running inside a Jupyter Notebook or Jupyter Lab."""
    try:
        get_ipython()  # pylint: disable=undefined-variable
        return True
    except Exception:  # pylint: disable=broad-exception-caught
        return False


def remove_lines(count=1):
    """Remove the last printed line(s) from the CLI."""
    if is_notebook_mode():
        # Jupyter
        # In Jupyter you can't clear a single line, only the entire cell output.
        clear_output(wait=True)
    else:
        # CLI
        while count > 0:
            count -= 1
            sys.stdout.write("\033[F")  # Move the cursor up one line
            sys.stdout.write("\033[K")  # Clear the line
            sys.stdout.flush()  # Flush the output buffer


def convertTuple(tup):
    # initialize an empty string
    if isinstance(tup, tuple):
        a_str = ""
        for item in tup:
            a_str = a_str + item
        return a_str

    return tup


def singular(string):
    return re.sub(r"s$", "", string)


def parse_path_tree(path_string):
    # Normalize the path string to use the appropriate separator for the current system.
    path = os.path.normpath(path_string)

    # Cut off first slash if it exists.
    if path[0] == "/":
        path = path[1:]

    # Remove file name and return the tree
    tree = path.split("/")[:-1]
    return tree


# Confirm promt for True or False Questions
def confirm_prompt(question: str = "", default=False) -> bool:
    reply = None
    while reply not in ("y", "n"):
        try:
            output_text(f"<yellow>{question}</yellow>", return_val=False)
            reply = input("(y/n): ").casefold()
            readline.remove_history_item(readline.get_current_history_length() - 1)
        except KeyboardInterrupt:
            print("\n")
            return default
    if reply == "y":
        return True
    return False


# Return boolean and formatted error message if other sessions exist.
def other_sessions_exist(cmd_pointer):
    from openad.app.global_var_lib import _meta_registry_session

    file_list = os.listdir(os.path.dirname(_meta_registry_session))
    try:
        file_list.remove("registry.pkl" + cmd_pointer.session_id)
    except Exception:  # pylint: disable=broad-exception-caught
        pass

    if len(file_list) > 0:
        output_error(msg("abort_clear_sessions"), return_val=False)
        return True
    else:
        return False


# Return user input.
def user_input(cmd_pointer, question):
    """
    Basically the same as input(), but with some extra styling and history disabled.
    """
    prompt = output_text(f"<yellow>{question}: </yellow>", return_val=True, jup_return_format="plain")
    text = input(prompt)
    return text


def user_secret(cmd_pointer, question):
    """
    Basically the same as getpass.getpass(), but with some extra styling and history disabled.
    """
    prompt = output_text(f"<yellow>{question}: </yellow>", return_val=True, jup_return_format="plain")
    text = getpass.getpass(prompt)
    return text


# Return list of available toolkit names.
def get_toolkits():
    folder_path = os.path.dirname(os.path.abspath(__file__)) + "/../user_toolkits"
    ignore_dirs = ["__pycache__", "readme"]
    toolkit_names = [
        name.upper()
        for name in os.listdir(folder_path)
        if os.path.isdir(os.path.join(folder_path, name)) and name not in ignore_dirs
    ]
    return toolkit_names


# Return boolean if toolkit is installed.
def is_toolkit_installed(toolkit_name, cmd_pointer=None):
    return cmd_pointer and toolkit_name and toolkit_name.upper() in cmd_pointer.settings["toolkits"]


# Validate a file path.
# - - -
# Add a default file extension when one is missing,
# or verify if the file extension is correct.
def validate_file_path(file_path: str, allowed_extensions: list, cmd_pointer):
    if not file_path:
        return

    default_extension = allowed_extensions[0]
    if len(file_path.split(".")) == 1:
        return file_path + "." + default_extension
    elif file_path.split(".")[-1].lower() not in allowed_extensions:
        output_error(msg("err_invalid_file_format", "csv"))
        return
    else:
        return file_path


# Ensure a file_path is kosher:
# - Make sure we won't override an existing file
# - Create folder structure if it doesn't exist yet
def ensure_file_path(file_path):
    if os.path.exists(file_path):
        # File already exists --> overwrite?
        if not confirm_prompt("The destination file already exists, overwrite?"):
            return False
    elif not os.path.isdir(os.path.dirname(file_path)):
        # Path doesn't exist --> create?
        if not confirm_prompt("The destination file path does not exist, create it?"):
            return False
        os.makedirs(os.path.dirname(file_path))
    return True


# Check is a port is open.
# Only used by next_avail_port below.
def is_port_open(host, port):
    import socket

    try:
        # Create a socket object and attempt to connect
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            s.settimeout(1)  # Set a timeout for the connection attempt
            s.connect((host, port))
        return False  # Port is occupied
    except (ConnectionRefusedError, socket.timeout):
        return True  # Port is available


# Return the next available port starting with 8024.
# This is used by the flask app launcher, we want to
# avoid a situation where multiple apps are trying to
# run on the same port.
def next_avail_port(port=8024, host="127.0.0.1"):
    # Not 127.0.0.1
    # - - -
    # In the context of interface binding, the address 127.0. 0.1
    # means that the server only listens to the loopback interface.
    # On the other hand, binding our server to the 0.0. 0.0 interface
    # means we want to accept traffic from all of the available interfaces.

    while not is_port_open(host, port):
        # print("port unavailable: ", port)
        port += 1
    return host, port


# Load python module from a dynamic path
def load_module_from_path(module_name, file_path):
    import importlib.util
    import sys

    try:
        spec = importlib.util.spec_from_file_location(module_name, file_path)
        module = importlib.util.module_from_spec(spec)
        sys.modules[module_name] = module
        spec.loader.exec_module(module)
        return module
    except Exception as err:
        # Silent fail - only enable this for debugging
        # output_error(f"load_module_from_path('{module_name}', {file_path})\n<soft>{err}</soft>")
        return None


# Print a terminal-wide separator
def print_separator(style=None, width=None, return_val=False):
    from openad.app.global_var_lib import GLOBAL_SETTINGS

    if GLOBAL_SETTINGS["display"] == "terminal" or GLOBAL_SETTINGS["display"] == None:
        cli_width = get_print_width(full=True)
        width = cli_width if not width or cli_width < width else width
        if style:
            return output_text(f"<{style}>{'-' * width}</{style}>", nowrap=True, return_val=return_val)
        else:
            return output_text(f"{'-' * width}", nowrap=True, return_val=return_val)


# Load a module or a module's function dynamically from a toolkit folder.
# This is a non-repo alt to `from foo import bar`
def load_tk_module(cmd_pointer, toolkit_name, lib_name, func_name=None):
    import importlib.util as ilu

    folder = cmd_pointer.toolkit_dir + f"/{toolkit_name}/{lib_name}.py"
    spec = ilu.spec_from_file_location(lib_name, folder)
    module = ilu.module_from_spec(spec)
    spec.loader.exec_module(module)
    if func_name:
        return getattr(module, func_name)
    else:
        return module


# Python equivalent of JavaScript's encodeURIComponent
def encode_uri_component(string):
    from urllib.parse import quote

    return quote(string.encode("utf-8"), safe="~()*!.'")


# Prettify a timestamp
def pretty_date(timestamp=None, style="log", include_time=True):
    # If no timestamp provided, use the current time
    if not timestamp:
        timestamp = time.time()

    # Choose the output format
    fmt = None
    if style == "log":
        fmt = "%d-%m-%Y"  # 07-01-2024
        if include_time:
            fmt += ", %H:%M:%S"  # 07-01-2024, 15:12:45
    elif style == "pretty":
        fmt = "%b %d, %Y"  # Jan 7, 2024
        if include_time:
            fmt += " at %H:%M"  # Jan 7, 2024 at 15:12
    else:
        output_error("Invalid style for pretty_date()")

    # Parse date/time string
    date_time = datetime.fromtimestamp(timestamp)
    return date_time.strftime(fmt)


# Prettify a number
def pretty_nr(nr, imperial=True):
    """
    Add commas to large numbers.

    Parameters
    ----------
    nr: int or float
        The number to format.
    imperial: bool
        Whether to use imperial formatting (commas) or not (spaces).

    Returns
    -------
    str:
        The formatted number as a string.
    """

    if nr is None and nr != 0:
        return None

    nr_split = str(nr).split(".")
    integer_str = nr_split[0]
    decimal_str = nr_split[1] if len(nr_split) > 1 else ""

    char = "," if imperial else " "
    output = re.sub(r"\B(?=(\d{3})+(?!\d))", char, integer_str)

    return output + (f".{decimal_str}" if decimal_str else "")


# Check if a variable (string or number) is numeric.
def is_numeric(str_or_nr):
    try:
        float(str_or_nr)
        return True
    except ValueError:
        return False


# Merge two lists of dictionaries while avoiding duplicates.
def merge_dict_lists(list1, list2):
    # Convert dictionaries to tuples of sorted items
    list1_tuples = [tuple(sorted(d.items())) for d in list1]
    list2_tuples = [tuple(sorted(d.items())) for d in list2]

    # Perform set operations to merge the lists while avoiding duplicates
    merged_tuples = list1_tuples + list(set(list2_tuples) - set(list1_tuples))

    # Convert the tuples back to dictionaries
    merged_list = [dict(t) for t in merged_tuples]

    return merged_list


# Get the available print width of the terminal.
def get_print_width(full=False):
    from openad.app.global_var_lib import GLOBAL_SETTINGS

    # Note: "api" can be removed after display() is removed from the %openadd magic command.

    # Notebook - fixed value
    if GLOBAL_SETTINGS["display"] == "notebook" or GLOBAL_SETTINGS["display"] == "api":
        return 120
    else:
        try:
            # Terminal full width
            if full:
                return shutil.get_terminal_size().columns

            # Terminal regular print width
            else:
                # We return the terminal width -10 so there's always room for
                # output with edge (5 chars) and some padding on the right.
                return min(shutil.get_terminal_size().columns - 10, GLOBAL_SETTINGS["max_print_width"])
        except Exception:  # pylint: disable=broad-exception-caught
            return GLOBAL_SETTINGS["max_print_width"]


# Style a boolean value in red or green
def style_bool(value):
    return (
        style(f"<success>{value}</success>")
        if value is True
        else style(f"<error>{value}</error>")
        if value is False
        else value
    )


def get_case_insensitive_key(dictionary, key_lowercase):
    """
    Get the key from a dictionary in a case-insensitive way.

    Parameters
    ----------
    dictionary: dict
        The dictionary to search in
    key_lowercase: str
        The key to search for

    Returns
    -------
    str:
        The matched key
    object:
        The value of the matched key
    """
    for key in dictionary:
        if key.lower() == key_lowercase.lower():
            return key, dictionary.get(key)

    return None, None


#
#
#
#
#
#
#
#
# NOT USED
#
#
#
#
#
#
#
#

# NOT USED
# Operating system


def get_platform():
    platforms = {"linux1": "Linux", "linux2": "Linux", "darwin": "OS X", "win32": "Windows"}
    if sys.platform not in platforms:
        return sys.platform
    return platforms[sys.platform]


# NOT USED
# Alternative to spinner that exits more gracefully on ctrl+c.
# Abandoning this for now because it's blocking the process.
def loader(text="", anim=["◐", "◓", "◑", "◒"], no_format=False, exit_msg=None, on_abort=None):
    import asyncio

    interval = 0.1
    sys.stdout.write("\033[?25l")  # Hide cursor

    async def _loop(text="", i=0, line_length=0):
        # print('>', line_length)
        a_str = anim[i]
        sys.stdout.write("\r" + " " * line_length + "\r")
        sys.stdout.flush()
        # Make text soft gray.
        text_formatted = text if no_format or not text else f"\u001b[90m{text}...\u001b[0m"
        line = "\r" + a_str + " " + text_formatted
        sys.stdout.write(line)  # This should return line_length but it doesn't for some reason.
        line_length = len(line)
        i = (i + 1) % len(anim)
        # time.sleep(interval) # Sync
        # _loop(text, i, line_length) # Sync
        await asyncio.sleep(interval)
        # await asyncio.run(_loop(text, i, line_length))

        try:
            await _loop(text, i, line_length)
            # _loop(text, i, line_length)
        except KeyboardInterrupt:
            print("\b \b \b" + "\b \b" * line_length + "\r")  # Clear the line.
            sys.stdout.write("\033[F")  # Move the cursor up one line
            if exit_msg:
                print(exit_msg)  # Move the cursor up one line
            if on_abort:
                on_abort()

    # _loop(text=text) # Sync
    asyncio.run(_loop(text=text))
    sys.stdout.write("\033[?25h")  # Show cursor


# NOT USED
# Clear the current line.
# Couldn't get this to work, because I can't clear the buffer.
# As a result, (eg.) when you try ctrl-c twice with n answer,
# the second time it will have the first time's response in the buffer,
# causing it to mess up the layout.
def clear_current_line():
    buffer = readline.get_line_buffer()
    line_length = len(buffer)
    readline.clear_history()
    eraser = "\b \b" * line_length
    sys.stdout.write(eraser)
    # readline.insert_text(' ')
    # print(len(buffer), buffer)
