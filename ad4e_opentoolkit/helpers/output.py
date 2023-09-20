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
from ad4e_opentoolkit.helpers.output_msgs import messages

# Importing our own plugins.
# This is temporary until every plugin is available as a public pypi package.
from ad4e_opentoolkit.plugins.style_parser import style, print_s, strip_tags, tags_to_markdown


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
    from ad4e_opentoolkit.helpers.general import is_notebook_mode
    notebook_mode = cmd_pointer.notebook_mode if cmd_pointer else is_notebook_mode()
    api_mode = cmd_pointer.api_mode if cmd_pointer else False
    return_val = notebook_mode if return_val is None else return_val

    if api_mode:
        # API

        return strip_tags(text)
    elif notebook_mode:

        # Jupyter
        if return_val:
            if jup_return_format == 'plain':
                return strip_tags(text)
            if jup_return_format == 'markdown_data':
                return Markdown(tags_to_markdown(text)).data
            else:
                return Markdown(tags_to_markdown(text))
        else:
            display(Markdown(tags_to_markdown(text)))
    else:
        # CLI
        if return_val:
            return style(text, **kwargs)
        else:
            print_s(text, **kwargs)


def _output_status(msg, status, cmd_pointer=None, return_val=None, pad_top=False, pad_btm=False, **kwargs):
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
    if isinstance(msg, tuple):
        msg1 = msg[0]
        msg2 = msg[1] if len(msg) == 2 else None
    else:
        msg1 = msg
        msg2 = None

    # Format secondary message.
    msg2 = f'\n<soft>{msg2}</soft>' if msg2 else ''

    # Set padding.
    pad = 0 if pad_top or pad_btm else 1

    # Print.
    return output_text(
        f'<{status}>{msg1}</{status}>{msg2}',
        cmd_pointer=cmd_pointer,
        return_val=return_val,
        pad=pad,
        pad_btm=pad_btm,
        pad_top=pad_top,
        **kwargs
    )


# Print or return error messages.
def output_error(msg, *args, **kwargs):
    """
    Wrapper around output_status() for error messages.
    """
    return _output_status(msg, 'error', *args, **kwargs)


# Print or return warning messages.
def output_warning(msg, *args, **kwargs):
    """
    Wrapper around output_status() for warning messages.
    """
    return _output_status(msg, 'warning', *args, **kwargs)


# Print or return success messages.
def output_success(msg, *args, **kwargs):
    """
    Wrapper around output_status() for success messages.
    """
    return _output_status(msg, 'success', *args, **kwargs)


# Print or return a table.
def output_table(data, cmd_pointer=None, headers=None, note=None, tablefmt='simple'):
    """
    Display a table:
    - CLI:      Print using tabulate with some custom home-made styling
    - Jupyter:  Return a pandas DataFrame
    - API:      Return a pandas DataFrame

    Parameters
    ----------
    data (list, required)
        An array of tuples, where each tuple is a row in the table.
    cmd_pointer (object):
        The run_cmd class instance which is used to determine the display context.
    headers (list):
        A list of strings, where each string is a column header.
    tablefmt (str):
        The table format used for tabulate (CLI only) - See https://github.com/astanin/python-tabulate#table-format
    """
    import shutil
    import pandas
    from tabulate import tabulate
    from ad4e_opentoolkit.helpers.general import is_notebook_mode
    notebook_mode = cmd_pointer.notebook_mode if cmd_pointer else is_notebook_mode()
    headers = [] if headers is None else headers
    is_df = isinstance(data, pandas.DataFrame)
    cli_width = shutil.get_terminal_size().columns

    # Turn potential tuples into lists.
    data = data if is_df else [list(row) for row in data]

    # Ensure the headers list matches the number of columns.
    col_count = data.shape[1] if is_df else len(data[0])

    if headers and len(headers) != col_count:
        output_warning(msg('table_headers_dont_match_columns', headers, col_count, split=True), cmd_pointer, return_val=False)
        headers = headers[:col_count] + ['(?)'] * max(0, col_count - len(headers))

    # - -
    # Return data in Jupyter.
    if notebook_mode is True:
        pandas.set_option('display.max_colwidth', None)
        # pandas.options.display.max_colwidth = 5
        # pandas.set_option('display.max_colwidth', 5)
        if (is_df):

            return data
        else:
            # Remove styling tags from headers.
            headers = list(map(lambda text: strip_tags(text), headers))

            # Remove styling tags from content.
            for i, row in enumerate(data):
                for j, cell in enumerate(row):
                    data[i][j] = strip_tags(cell)

            return pandas.DataFrame(data, columns=headers)

    # - -
    # Display data in terminal.
    if (is_df):
        table = tabulate(data, headers="keys", tablefmt=tablefmt, showindex=False, numalign="left")
    else:
        # Parse styling tags.
        for i, row in enumerate(data):
            for j, cell in enumerate(row):
                data[i][j] = style(cell, nowrap=True)

        table = tabulate(data, headers=headers, tablefmt=tablefmt, showindex=False, numalign="left")

    # Crop table if it's wider than the terminal.
    max_row_length = max(list(map(lambda row: len(row), table.splitlines())))
    if max_row_length > cli_width:
        for i, line in enumerate(table.splitlines()):
            if i == 1:
                table = table.replace(line, line[:cli_width])
            elif len(line) > cli_width:
                # updated with reset \u001b[0m for color tags which may be found later
                table = table.replace(line, line[:cli_width - 3] + '...\u001b[0m')

    # Make line yellow
    lines = table.splitlines()
    lines[1] = style(f'<yellow>{lines[1]}</yellow>', nowrap=True)

    # List next commands
    msg = (
        '<cmd>show full table</cmd>',
        '<cmd>remove table columns</cmd>',
        '<cmd>edit table</cmd>',
        '<cmd>save table as \'<filename.csv>\'</cmd>',
    )
    lines.append('\n<soft>Next up you can: </soft>' + ' / '.join(msg) + ' <green># Coming soon</green>')

    # Add footnote.
    if note:
        lines.append(f'\n<soft>{note}</soft>')
    table = '\n'.join(lines)

    # Print.
    print_s(table, pad=2, nowrap=True)


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
            msg_name = '\n'.join(msg_name)
    return msg_name
