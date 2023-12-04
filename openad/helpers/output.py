"""
This library lets us print styled content in the CLI and Jupyter.
Use this instead of print() which should always be avoided.

Under the Hood
--------------
- The style_parser plugin takes care of parsing xml tags into ANSI
  escape codes, which lets us print different colors in the CLI.
  It also takes care of some extra layout fluff, like padding,
  indentation, etc. The main functions are style() and print_s()
  which return & print styled text respectively. Please refer to
  the style_parser documentation for more details on how to style text.
- output_text() takes care of checking what environment we're in,
  (CLI, Jupyter or API) and will either print or return the text
  styled according to the relevant display context. In the CLI it
  will use the style_parser, in Jupyter it will convert the XML to
  Markdown, and in the API it will strip the tags and return plain text.

Usage
-----
output_text()
    This is the main output function and replaces print().
    Simple usage: output_text('Hello world', cmd_pointer, pad=2)
    Note: you can pass additional parameters which will be passed
    onto the style_parser, eg. pad=2. See documentation for style().
output_error()
output_warning()
output_success()
    These are all wrappers around _output_status().
_output_status()
    This is a wrapper around output_text() which takes care of
    some templated styling to make sure error/waring/success outputs
    are always treated consistently.
output_table()
    Displays a table in the CLI, or returns a panda dataframe when
    called from Jupyter or the API.

Examples:
---------
output_error(msg('err_login', toolkit_name, err, split=True), cmd_pointer)
input(output_text('<yellow>Hello</yellow>', jup_return_format='plain')
output_text(msg('workspace_description', workspace_name, description), cmd_pointer, pad=1, edge=True)
"""

from IPython.display import Markdown, display
from openad.helpers.output_msgs import messages


# NOTE (moe) - I'm not able to do import memory from the core director
# Meanwhile this is obsolete, but I still don't have an explanation.
# - - -
# from openad.app.memory import memory # This works
# from openad.core.memory import memory # This crashes the app
# import openad.core.memory as memory # This crashes the app
# from openad.core.test import foo # This also crashes the app, when test.py has `foo = 123`


# Importing our own plugins.
# This is temporary until every plugin is available as a public pypi package.
from openad.plugins.style_parser import style, print_s, strip_tags, tags_to_markdown


# Print or return styled text.
def output_text(text, cmd_pointer=None, return_val=None, jup_return_format=None, **kwargs):
    """
    Print or return styled text according to the relevant display context (API/Jupyter/CLI).

    Parameters:
    -----------
    text (str, required):
        The text to display. This can be plain text, but in most cases we want to
        store the message in output_msgs.py and request it using msg('<msg_name>').
    cmd_pointer (object):
        The run_cmd class instance which is used to determine the display context.
        While it is technically optional, it is recommended to always pass it.
        We can still figure out the CLI/Jupyter display context without it, but
        we won't be able to do that for the API.
    return_val (None/bool):
        Return the text instead of displaying it.
        The default (None) results in True for Jupyter, False for CLI.
    jup_return_format (None/'plain'/'markdown_data'):
        By default, we return a markdown object in Jupyter, but sometimes
        we just need plain text (eg. for input() or the Halo spinner).
        When we need to append the markdown content to another string
        pre-formatting, we want to return the markdown data instead.
        The latter is the case with the plash.py installation notice.
    """
    from openad.helpers.general import is_notebook_mode

    notebook_mode = cmd_pointer.notebook_mode if cmd_pointer else is_notebook_mode()
    api_mode = cmd_pointer.api_mode if cmd_pointer else False
    return_val = notebook_mode if return_val is None else return_val

    # API
    if api_mode:
        return strip_tags(text)

    # Jupyter
    elif notebook_mode:
        if return_val:
            if jup_return_format == "plain":
                return strip_tags(text)
            if jup_return_format == "markdown_data":
                return Markdown(tags_to_markdown(text)).data
            else:
                return Markdown(tags_to_markdown(text))
        else:
            display(Markdown(tags_to_markdown(text)))

    # CLI
    else:
        if return_val:
            return style(text, **kwargs)
        else:
            print_s(text, **kwargs)


def _output_status(message, status, cmd_pointer=None, return_val=None, pad=1, pad_top=False, pad_btm=False, **kwargs):
    """
    Print or return styled error/warning/success message according
    to the relevant display context (API/Jupyter/CLI).

    This is simply a style wrapper around output_text().
    By default text is displayed with pad=1, unless either
    pad_btm/pad_top is set.

    Parameters
    ----------
    msg (str/tuple, required):
        The message to display. If a tuple is passed, the first item
        is the main message in red/orange/green, and the second is a
        secondary message in soft grey.
    status (str, required):
        The status type. One of 'error', 'warning' or 'success'.
    cmd_pointer (object):
        The run_cmd class instance which is used to determine the display context.
    return_val (None/bool):
        Return the styled text instead of printing it. This will override
        the default which is True for Jupyter, False for CLI.
    pad_top (bool):
        Add padding to the top of the message only, instead of the default
        of 1 line of padding both top and bottom (see style() documentation).
    pad_btm (bool):
        Add padding to the bottom of the message only, instead of the default
        of 1 line of padding both top and bottom (see style() documentation).
    kwargs:
        Additional parameters to pass onto the style_parser.
    """

    # Check if there's a secondary message.
    if isinstance(message, tuple):
        msg1 = message[0]
        msg2 = message[1] if len(message) == 2 else None
    else:
        msg1 = message
        msg2 = None

    # Format secondary message.
    msg2 = f"\n<soft>{msg2}</soft>" if msg2 else ""

    # Set padding.
    pad = 0 if pad_top or pad_btm else pad

    # Print.
    return output_text(
        f"<{status}>{msg1}</{status}>{msg2}",
        cmd_pointer=cmd_pointer,
        return_val=return_val,
        pad=pad,
        pad_btm=pad_btm,
        pad_top=pad_top,
        **kwargs,
    )


# Print or return error messages.
def output_error(message, *args, **kwargs):
    """
    Wrapper around output_status() for error messages.
    """
    return _output_status(message, "error", *args, **kwargs)


# Print or return warning messages.
def output_warning(message, *args, **kwargs):
    """
    Wrapper around output_status() for warning messages.
    """
    return _output_status(message, "warning", *args, **kwargs)


# Print or return success messages.
def output_success(message, *args, **kwargs):
    """
    Wrapper around output_status() for success messages.
    """
    return _output_status(message, "success", *args, **kwargs)


# Print or return a table.
def output_table(table, cmd_pointer=None, is_data=False, headers=None, note=None, tablefmt="simple"):
    """
    Display a table:
    - CLI:      Print using tabulate with some custom home-made styling
    - Jupyter:  Return a pandas DataFrame
    - API:      Return a pandas DataFrame

    Parameters
    ----------
    data (list, required)
        A dataframe or an array of tuples, where each tuple is a row in the table.
    cmd_pointer (object):
        The run_cmd class instance which is used to determine the display context.
    is_data (bool):
        This enables the follow-up commands to open/edit/save the table data.
        Some tables are just displaying information and don't need this (eg workspace list.)
    headers (list):
        A list of strings, where each string is a column header.
    note (str):
        A footnote to display at the bottom of the table.
    tablefmt (str):
        The table format used for tabulate (CLI only) - See https://github.com/astanin/python-tabulate#table-format
    """
    import shutil
    import pandas
    from tabulate import tabulate
    from openad.helpers.general import is_notebook_mode

    notebook_mode = cmd_pointer.notebook_mode if cmd_pointer else is_notebook_mode()
    headers = [] if headers is None else headers
    is_df = isinstance(table, pandas.DataFrame)
    cli_width = shutil.get_terminal_size().columns

    # Abort when table is empty.
    if _is_empty_table(table, is_df):
        output_warning(msg("table_is_empty"), cmd_pointer, return_val=False)
        return

    # Turn potential tuples into lists.
    table = table if is_df else [list(row) for row in table]

    # Ensure the headers list matches the number of columns.
    col_count = table.shape[1] if is_df else len(table[0])

    if headers and len(headers) != col_count:
        output_warning(
            msg("table_headers_dont_match_columns", headers, col_count, split=True), cmd_pointer, return_val=False
        )
        headers = headers[:col_count] + ["(?)"] * max(0, col_count - len(headers))

    # Enable follow-up commands.
    if is_data:
        if not cmd_pointer:
            raise Exception("cmd_pointer is required in display_data() to enable follow-up commands.")

        # Store data in memory so we can access it with follow-up commands.
        cmd_pointer.memory.store(table)

    # - -
    # Format data for Jupyter.
    if notebook_mode is True:
        pandas.set_option("display.max_colwidth", None)
        # pandas.options.display.max_colwidth = 5
        # pandas.set_option('display.max_colwidth', 5)
        if is_df:
            # return data %%
            pass
        else:
            # Remove styling tags from headers.
            if headers:
                headers = list(map(lambda text: strip_tags(text), headers))

            # Remove styling tags from content.
            for i, row in enumerate(table):
                for j, cell in enumerate(row):
                    table[i][j] = strip_tags(cell)

            # return pandas.DataFrame(data, columns=headers) %%
            table = pandas.DataFrame(table, columns=headers)

    # - -
    # Format data for terminal.
    else:
        if is_df:
            table = tabulate(table, headers="keys", tablefmt=tablefmt, showindex=False, numalign="left")
        else:
            # Parse styling tags.
            for i, row in enumerate(table):
                for j, cell in enumerate(row):
                    table[i][j] = style(cell, nowrap=True)

            table = tabulate(table, headers=headers, tablefmt=tablefmt, showindex=False, numalign="left")

        # Crop table if it's wider than the terminal.
        max_row_length = max(list(map(lambda row: len(row), table.splitlines())))
        if max_row_length > cli_width:
            for i, line in enumerate(table.splitlines()):
                if i == 1:
                    table = table.replace(line, line[:cli_width])
                elif len(line) > cli_width:
                    # updated with reset \u001b[0m for color tags which may be found later
                    table = table.replace(line, line[: cli_width - 3] + "...\u001b[0m")

        # Make line yellow.
        table = table.splitlines()
        table[1] = style(f"<yellow>{table[1]}</yellow>", nowrap=True)
        table = "\n".join(table)

    # Display footnote.
    footnote = ""

    # --> Optional follow-up commands.
    if is_data:
        message = (
            "<cmd>result open</cmd>",
            "<cmd>edit</cmd>",
            "<cmd>copy</cmd>",
            "<cmd>display</cmd>",
            "<cmd>save [as '<filename.csv>']</cmd>",
        )
        footnote += "<soft>Next up, you can run: </soft>" + "/".join(message)

    # --> Optional custom note.
    if note:
        footnote += f"\n<soft>{note}</soft>"

    # Output
    if notebook_mode is True:
        if footnote:
            output_text(footnote, cmd_pointer, return_val=False)
        return table
    else:
        if footnote:
            output = table + "\n\n" + footnote
        else:
            output = table
        print_s(output, pad=2, nowrap=True)


# Check whether table data is empty.
def _is_empty_table(data, is_df):
    if is_df:
        return data.empty
    elif not data:
        return True
    else:
        is_empty = True
        for col in data:
            if col:
                is_empty = False
                break
        return is_empty


# Procure a display message from output_msgs.py.
def msg(msg_name, *args, split=False):
    """
    Fetches the correct output message from the messages dictionary.

    Output messages are stored in output_msgs.py in three different
    formats, or a combination of them:
    - String: for simple static messages
    - Lambda function: for messages with variables
    - Tuple: for messages with multiple lines

    Parameters
    ----------
    msg_name (str):
        The name of the message to fetch.
    args:
        Any number of variables that are required for the lambda function.
    split (bool):
        Instead of parsing a tuple as different lines of the same string,
        you can return it as a list of separate messages. This is used for
        the error/warning/success messages, where the first message is the
        main message, and the second message is a secondary message.
    """
    msg_name = messages[msg_name]
    if callable(msg_name):
        msg_name = msg_name(*args)
    if isinstance(msg_name, tuple):
        if not split:
            # For output_error/warning/success we sometimes need to send
            # None as second message, eg. load_toolkit_description().
            msg_name = [x for x in msg_name if x is not None]
            msg_name = "\n".join(msg_name)
    return msg_name
