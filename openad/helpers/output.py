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
    Simple usage: output_text('Hello world', pad=2)
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
output_error(msg('err_login', toolkit_name, err))
input(output_text('<yellow>Hello</yellow>', jup_return_format='plain')
output_text(msg('workspace_description', workspace_name, description), pad=1, edge=True)
"""

import shutil
import pandas
from tabulate import tabulate
from IPython.display import Markdown
from openad.helpers.output_msgs import msg

# Importing our own plugins.
# This is temporary until every plugin is available as a public pypi package.
from openad.plugins.style_parser import style, print_s, strip_tags, tags_to_markdown


# NOTE (moe) - I'm not able to do import memory from the core director
# Meanwhile this is obsolete, but I still don't have an explanation.
# - - -
# from openad.app.memory import memory # This works
# from openad.core.memory import memory # This crashes the app
# import openad.core.memory as memory # This crashes the app
# from openad.core.test import foo # This also crashes the app, when test.py has `foo = 123`


# Print or return styled text.
def output_text(message, return_val=None, jup_return_format=None, **kwargs):
    """
    Print or return styled text according to the relevant display context (API/Jupyter/CLI).

    Parameters:
    -----------
    text (str, required):
        The text to display. This can be plain text, but in most cases we want to
        store the message in output_msgs.py and request it using msg('<msg_name>').
    return_val (None/bool):
        Returns the text instead of displaying it.
        The default (None) results in True for Jupyter, False for CLI.
        NOTE: This is confusing and will be refactored so return is always consistent regardless of environment.
    jup_return_format (None/'plain'/'markdown_data'):
        By default, we return a markdown object in Jupyter, but sometimes
        we just need plain text (eg. for input() or the Halo spinner).
        When we need to append the markdown content to another string
        pre-formatting, we want to return the markdown data instead.
        The latter is the case with the plash.py installation notice.
    kwargs:
        Additional parameters for style_parser.
    """

    # Imported here to avoid circular imports.
    from openad.app.global_var_lib import GLOBAL_SETTINGS

    DISPLAY = GLOBAL_SETTINGS["display"]
    return_val = DISPLAY == "notebook" if return_val is None else return_val

    # When the message is a list of strings, the first string
    # will be printed regularly and subsequent strings will be
    # printed soft gray.
    # - - -
    # This is not really used but it's to ensure output
    # is consistent between output_text and output_error etc.
    if isinstance(message, list):
        message = "\n".join([f"<soft>{string}</soft>" if i > 0 else string for i, string in enumerate(message)])

    # API
    if DISPLAY == "api":
        return strip_tags(message)

    # Jupyter
    elif DISPLAY == "notebook":
        if return_val:
            if jup_return_format == "plain":
                return strip_tags(message)
            if jup_return_format == "markdown_data":
                return Markdown(tags_to_markdown(message)).data
            else:
                return Markdown(tags_to_markdown(message))
        else:
            DISPLAY(Markdown(tags_to_markdown(message)))

    # CLI
    elif DISPLAY == "terminal" or DISPLAY == None:
        if return_val:
            return style(message, **kwargs)
        else:
            print_s(message, **kwargs)


def _output_status(message, status, pad=1, pad_top=None, pad_btm=None, *args, **kwargs):
    """
    Assure consistent styling for error/warning/success messages.

    This is simply a style wrapper around output_text() which takes care of:
    - Displaying message in correct color
    - Add 1 line of spacing before and after, unless specified otherwise

    Parameters
    ----------
    status (None/'error'/'warning'/'success'):
        Controls what color (red/yellow/green) the status message is assigned.
    pad (int, default=1):
        Used by style_parser: Adds padding to the top and bottom of the message.
        Gets ignored when either pad_top or pad_btm is set.
    pad_top (bool, default=0):
        Used by style_parser: Ignores the pad value and adds padding to the top of the message only.
    pad_btm (bool, default=0):
        Used by style_parser: Ignores the pad value and adds padding to the bottom of the message only.
    """

    # When the message is a list of strings, the first string
    # will be printed in the status color (red/yellow/green)
    # and subsequent strings will be printed soft gray.
    if isinstance(message, list):
        message = "\n".join(
            [
                f"<soft>{string}</soft>" if i > 0 else f"<{status}>{string}</{status}>"
                for i, string in enumerate(message)
            ]
        )
    else:
        message = f"<{status}>{message}</{status}>"

    # Set default padding of 1, unless pad_top or pad_btm is set.
    # These parameters are consumed by style_parser, hence they are
    # wrapped into kwargs to keep output_text simple.
    pad = 0 if type(pad_top) == int or type(pad_btm) == int else pad
    kwargs["pad"] = pad
    kwargs["pad_btm"] = pad_btm
    kwargs["pad_top"] = pad_top

    return output_text(message, *args, **kwargs)


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
def output_table(table, is_data=True, headers=None, note=None, tablefmt="simple"):
    """
    Display a table:
    - CLI:      Print using tabulate with some custom home-made styling
    - Jupyter:  Return a pandas DataFrame
    - API:      Return a pandas DataFrame

    Parameters
    ----------
    table (dataframe/list, required)
        A dataframe or a list of tuples, where each tuple is a row in the table.
    is_data (bool):
        Enables the follow-up commands (open/edit/save/copy). Some tables are
        just displaying information and don't need them (eg workspace list.)
    headers (list):
        A list of strings, where each string is a column header.
    note (str):
        A footnote to display at the bottom of the table.
    tablefmt (str):
        The table format used for tabulate (CLI only) - See https://github.com/astanin/python-tabulate#table-format
    """

    # Imported here to avoid circular imports.
    from openad.app.global_var_lib import MEMORY
    from openad.app.global_var_lib import GLOBAL_SETTINGS

    DISPLAY = GLOBAL_SETTINGS["display"]

    headers = [] if headers is None else headers
    is_df = isinstance(table, pandas.DataFrame)
    cli_width = shutil.get_terminal_size().columns

    # Abort when table is empty.
    if _is_empty_table(table, is_df):
        output_warning(msg("table_is_empty"), return_val=False)
        return

    # Turn potential tuples into lists.
    table = table if is_df else [list(row) for row in table]

    # Ensure the headers list matches the number of columns.
    col_count = table.shape[1] if is_df else len(table[0])

    if headers and len(headers) != col_count:
        output_warning(msg("table_headers_dont_match_columns", headers, col_count), return_val=False)
        headers = headers[:col_count] + ["(?)"] * max(0, col_count - len(headers))

    # Enable follow-up commands.
    if is_data:
        MEMORY.store(table)

    # - -
    # Format data for Jupyter.
    if DISPLAY == "notebook":
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
            else:
                headers = None

            # Remove styling tags from content.
            for i, row in enumerate(table):
                for j, cell in enumerate(row):
                    table[i][j] = strip_tags(cell)

            table = pandas.DataFrame(table, columns=headers)

    # - -
    # Format data for terminal.
    elif DISPLAY == "terminal" or DISPLAY == None:
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

        # Make header line yellow.
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
        footnote += f"<soft>{note}</soft>"

    # Output
    if DISPLAY == "notebook":
        if footnote:
            output_text(footnote, return_val=False)
        return table
    elif DISPLAY == "terminal" or DISPLAY == None:
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
