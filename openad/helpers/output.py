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

import math
import pandas
from tabulate import tabulate
from IPython.display import Markdown, display
from openad.helpers.output_msgs import msg
import openad.helpers.general as helpers_general

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

    return_val = GLOBAL_SETTINGS["display"] == "notebook" if return_val is None else return_val

    # When the message is a list of strings, the first string
    # will be printed regularly and subsequent strings will be
    # printed soft gray.
    # - - -
    # This is not really used but it's to ensure output
    # is consistent between output_text and output_error etc.
    if isinstance(message, list):
        message = "\n".join([f"<soft>{string}</soft>" if i > 0 else string for i, string in enumerate(message)])

    # API
    if GLOBAL_SETTINGS["display"] == "api":
        return strip_tags(message)

    # Jupyter
    elif GLOBAL_SETTINGS["display"] == "notebook":
        if return_val:
            if jup_return_format == "plain":
                return strip_tags(message)
            if jup_return_format == "markdown_data":
                return Markdown(tags_to_markdown(message)).data
            else:
                return Markdown(tags_to_markdown(message))
        else:
            display(Markdown(tags_to_markdown(message)))

    # CLI
    elif GLOBAL_SETTINGS["display"] == "terminal" or GLOBAL_SETTINGS["display"] == None:
        if return_val:
            return style(message, **kwargs)
        else:
            print_s(message, **kwargs)


def _output_status(message, status, pad=None, pad_top=None, pad_btm=None, *args, **kwargs):
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
def output_table(
    table,
    is_data=True,
    headers=None,
    show_index=False,
    note=None,
    return_val=None,
    pad=1,
    pad_btm=None,
    pad_top=None,
    _display=None,
    **kwargs,  # For tabulate options
):
    """
    Display a table:
    - CLI:      Print using tabulate with some custom home-made styling
    - Jupyter:  Return a pandas DataFrame
    - API:      Return a pandas DataFrame

    Parameters
    ----------
    table (dataframe/list, required)
        A dataframe or a list of lists/tuples, where each list/tuple is a row
        in the table.
    is_data (bool, default=True):
        Enables the follow-up commands (open/edit/save/copy). Some tables are
        just displaying information and don't need them (eg workspace list.)
    headers (list):
        A list of strings, where each string is a column header.
        Ignored when table is a dataframe.
    show_index (bool, default=False):
        Display the index column in Jupyter.
    note (str):
        A footnote to display at the bottom of the table.
    return_val (None/bool):
        Returns the table instead of displaying it.
        The default (None) results in True for Jupyter and API, False for CLI.
        NOTE: This is confusing and will be refactored so return is always
        consistent regardless of environment.
    pad (int, default=1):
    pad_top (int):
    pad_btm (int):
        Mimics the behavior of the style parser consistent with output_text etc.
        These parameters control how many empty lines are displayed before and
        after the table. If pad_btm or pad_top is set, pad will be ignored.
    _display (bool):
        Used for testing only. When running output_table outside of the OpenAD
        context, directly from a Python file or a Notebook, GLOBAL_SETTINGS["display"]
        will always be none, so this way we can set it manually.
    """

    # Imported here to avoid circular imports.
    from openad.app.global_var_lib import MEMORY
    from openad.app.global_var_lib import GLOBAL_SETTINGS

    # A dataframe can be styled, eg. df.style.set_properties(**{"text-align": "left"})
    # but a styled dataframe comes in as a styler object, which breaks functionality.
    # So we extract the dataframe from the styler object, and then override the styler
    # object's original data with the manipulated data before printing.
    is_df_styler_obj = hasattr(pandas.io.formats, "style") and isinstance(table, pandas.io.formats.style.Styler)
    if is_df_styler_obj:
        table_styler = table
        table = table.data

    is_df = isinstance(table, pandas.DataFrame)
    cli_width = helpers_general.get_print_width(full=True)

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

    # Prioritize output_table headers over dataframe columns.
    if is_df and headers:
        table.columns = headers

    # - -
    # Format data for Jupyter.
    if GLOBAL_SETTINGS["display"] in ["api", "notebook"] or _display == "notebook":
        pandas.set_option("display.max_colwidth", None)
        if is_df:
            pass
        else:
            # Remove styling tags from headers.
            if headers:
                headers = list(map(lambda text: strip_tags(text), headers))

            # Remove styling tags from content.
            for i, row in enumerate(table):
                for j, cell in enumerate(row):
                    table[i][j] = strip_tags(cell)

            table = pandas.DataFrame(table, columns=headers)

    # - -
    # Format data for terminal.
    elif GLOBAL_SETTINGS["display"] == "terminal" or _display == "terminal":
        if is_df:
            # By default, display the columns from the dataframe.
            tabulate_headers = "keys"

            # Remove the default numeric header if no columns are set.
            if not headers:
                if table.columns.values.tolist()[0] == 0 and table.columns.values.tolist()[1] == 1:
                    tabulate_headers = []
                else:
                    headers = table.columns.values.tolist()

            table = tabulate(
                table, headers=tabulate_headers, tablefmt="simple", showindex=False, numalign="right", **kwargs
            )
        else:
            # Parse styling tags.
            for i, row in enumerate(table):
                for j, cell in enumerate(row):
                    table[i][j] = style(cell, nowrap=True)

            headers = [] if headers is None else headers
            table = tabulate(table, headers=headers, tablefmt="simple", showindex=False, numalign="right", **kwargs)

        # Crop table if it's wider than the terminal.
        max_row_length = max(list(map(lambda row: len(row), table.splitlines())))
        if max_row_length > cli_width:
            for i, line in enumerate(table.splitlines()):
                if i == 1:
                    table = table.replace(line, line[:cli_width])
                elif len(line) > cli_width:
                    # updated with reset \u001b[0m for color tags which may be found later
                    table = table.replace(line, line[: cli_width - 3] + "...\u001b[0m")

    # Display footnote.
    footnote = ""

    # --> Optional follow-up commands.
    if is_data:
        message = (
            "<cmd>result open</cmd>",
            "<cmd>edit</cmd>",
            "<cmd>copy</cmd>",
            "<cmd>display</cmd>",
            "<cmd>as dataframe</cmd>",
            "<cmd>save [as '<filename.csv>']</cmd>",
        )
        footnote += "<soft>Next up, you can run: </soft>" + "/".join(message)

    # --> Optional custom note.
    if note:
        footnote += f"<soft>{note}</soft>"

    # Output for Jupyter
    if GLOBAL_SETTINGS["display"] in ["api", "notebook"] or _display == "notebook":
        # If the passed dataframe is a styler object, we extracted the
        # dataframe earlier, and now re-assign the manipulated dataframe
        # back to the styler object.
        if is_df_styler_obj:
            table_styler.data = table
            table = table_styler

            # Hide index column in Jupyter, unless we want it.
            if not show_index:
                table.hide(axis="index")
        else:
            if not show_index:
                table = table.style.hide(axis="index")

        if footnote:
            output_text(footnote, return_val=False)

        if GLOBAL_SETTINGS["display"] == "api":
            if isinstance(table, pandas.io.formats.style.Styler):
                return table.data
            else:
                return table
        elif return_val is True or return_val is None:
            return table
        else:
            display(table)

    # Output for terminal
    elif GLOBAL_SETTINGS["display"] == "terminal" or _display == "terminal":
        # Return table.
        if return_val is True:
            return table

        # Print table.
        else:
            if pad_top:
                print("\n".join([""] * pad_top))
            elif pad and not pad_btm:
                print("\n".join([""] * pad))

            _paginated_output(table, headers=headers, exit_msg=footnote)

            if pad_btm:
                print("\n".join([""] * pad_btm))
            elif pad and not pad_top:
                print("\n".join([""] * pad))


def _paginated_output(table, headers=None, exit_msg=None):
    from openad.helpers.general import remove_lines, print_separator
    import shutil

    output_lines = table.split("\n")

    # Calculate table width.
    table_width = 0
    for line in output_lines:
        if len(line) > table_width:
            table_width = len(line)

    # Divide output lines in 3 sections.
    lines_header = []
    lines_body = []

    # Isolate header lines.
    header_height = 0
    if headers:
        for column in headers:
            height = len(column.splitlines())
            if height > header_height:
                header_height = height
    lines_header = output_lines[: header_height + 1]

    # Isolate body lines.
    lines_body = output_lines[header_height + 1 :]

    # Color separator(s) yellow.
    sep = lines_header[len(lines_header) - 1]
    sep = style(f"<yellow>{sep}</yellow>", nowrap=True)
    lines_header = lines_header[:-1] + [sep]
    if not headers:
        sep2 = lines_body[len(lines_body) - 1]
        sep2 = style(f"<yellow>{sep2}</yellow>", nowrap=True)
        lines_body = lines_body[:-1] + [sep2]

    row_cursor = 0
    skip = 0
    while row_cursor < len(output_lines):
        # Minus 3 lines for bottom ellipsis, separator and pagination.
        max_rows = shutil.get_terminal_size().lines - len(lines_header) - 3

        # Print header + separator.
        if row_cursor == 0:
            row_cursor += len(lines_header)
        else:
            print("")
        print("\n".join(lines_header))

        # Print body.
        i = 0
        while i < max_rows and row_cursor < len(output_lines):
            print(lines_body[skip + i])
            i = i + 1
            row_cursor = row_cursor + 1
        skip = skip + i

        # Print pagination.
        if row_cursor < len(output_lines):
            total_pages = math.ceil(len(output_lines) / max_rows)
            page = math.floor(row_cursor / max_rows)

            try:
                print("...")
                print_separator("yellow", table_width)
                input(
                    style(
                        f"<yellow>Page {page}/{total_pages}</yellow> - Press <reverse> ENTER </reverse> to load the next page, or <reverse> ctrl+C </reverse> to exit."
                    )
                )

                # [ WAIT FOR INPUT ]

                remove_lines(3)
            except KeyboardInterrupt:
                print()  # Start a new line.
                remove_lines(2)
                if exit_msg:
                    print_s("\n" + exit_msg)
                return

        # Print exit message.
        else:
            if exit_msg:
                print_s("\n" + exit_msg)


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
