"Pretty print data in columns."

import math
from openad.helpers.general import get_print_width


def list_columns(
    values: list,
    print_width: int = None,
    col_width: int = 40,
    is_truncated: list = False,
    indent: int = 0,
):
    """
    Display a list of values in columns.
    """

    if not print_width:
        print_width = get_print_width()
    print_width = print_width - indent
    output = ""
    line_nr = 0

    # Remove blank values
    values = [f"{value}" for value in values if value]

    # The number of columns
    col_count = print_width // col_width

    # The length of each column
    col_length = math.ceil(len(values) / col_count)

    # Remove enough values so the total is dividable by column count
    values = values[: len(values) - (len(values) % col_count)]

    # Truncate the values that are longer than the column width.
    values = [f"{indent * ' '}{value[:col_width - 5]}..." if len(value) > col_width - 2 else value for value in values]

    # Pad each value so it fits the column width
    values = [f"{indent * ' '}{value:<{col_width}}" for value in values]

    # Compile the output
    for val in values:
        if line_nr > 0 and line_nr % col_count == 0:
            output = output + "\n"
        output = output + val
        line_nr = line_nr + 1 % col_length

    # Add elipsis to the end of each column if the output is truncated
    if is_truncated:
        output = output + "\n" + (f"{'<soft>...</soft>':<{col_width + 13}}" * col_count)

    return output


def key_val_columns(
    items_dict: dict,
    print_width: int = None,
    col_width: int = 40,
    ignore_keys: list = None,
    indent: int = 0,
):
    """
    Display a dictionary's keys + values in columns.
    """

    if not print_width:
        print_width = get_print_width()
    print_width = print_width - indent
    output = ""
    gap = 4

    # The number of columns
    col_count = print_width // col_width

    # Gap between the columns
    col_width = col_width - gap

    items = []
    for key, val in items_dict.items():
        # Ignore keys
        if ignore_keys and key in ignore_keys:
            continue

        key_len = len(str(key))
        val_len = len(str(val))

        # Gray out blank values and replace with '-'
        if val != 0 and not val:
            val = "<soft>-</soft>"
            val_len = 1
            key_color = "soft"
        else:
            key_color = "cyan"

        # Key + value on one line
        if key_len + 2 + val_len <= col_width:
            spacer = "<soft>" + ("." * (col_width - key_len - val_len - 1)) + "</soft>"
            key_val = f"{(indent * ' ')}<{key_color}>{key}:</{key_color}>{spacer}{val}"
            items.append(key_val + (" ") * gap)

        # Key + value split over two lines,
        # with possibly truncated value
        else:
            spacer_key = "<soft>" + ("." * (col_width - key_len - 1)) + "</soft>"
            spacer_val = "<soft>" + ("." * (col_width - val_len)) + "</soft>"
            val = f"{val[:col_width - 3]}..." if val_len > col_width else val
            key_val_1 = f"{(indent * ' ')}<{key_color}>{key}</{key_color}>:{spacer_key}"
            key_val_2 = f"{(indent * ' ')}{spacer_val}{val}"
            items.append(key_val_1 + (" ") * gap)
            items.append(key_val_2 + (" ") * gap)

    # The length of each column
    col_length = math.ceil(len(items) / col_count)

    # Arrange items into columns
    table = []
    for j, item in enumerate(items):
        if j % col_length == 0:
            table.append([])
        table[-1].append(item)

    # Reorder the items to the print top to bottom, left to right
    # instead of left to right top to bottom
    table = _transpose_table(table)

    # Compile the output
    for row in table:
        output = output + "".join(row) + "\n"

    return output


def _transpose_table(table):
    """
    Transpose a table (list of lists).

    Input:
    [
      [a,b,c,d,e],
      [f,g,h,i,j],
      [k,l,m]
    ]

    Output:
    [
      [a,f,k],
      [b,g,l],
      [c,h,m],
      [d,i,None],
      [e,j,None]
    ]
    """
    col_len = len(table[0])
    table[-1].extend([""] * (col_len - len(table[-1])))
    return list(zip(*table))


def key_val_full(
    items_dict: dict,
    print_width: int = None,
    ignore_keys: list = None,
    indent: int = 0,
):
    """
    Display a dictionary's keys + values covering the full width.
    """

    if not print_width:
        print_width = get_print_width()

    items = []
    for key, val in items_dict.items():
        # Ignore keys
        if ignore_keys and key in ignore_keys:
            continue

        val_len = len(str(val))

        # Gray out blank values and replace with '-'
        if val != 0 and not val:
            val = "<soft>-</soft>"
            val_len = 1
            key_color = "soft"
        else:
            key_color = "cyan"

        val = f"{val[:print_width - len(key) - 2 - 3]}..." if val_len > print_width else val
        key_val = f"{(indent * ' ')}<{key_color}>{key}:</{key_color}> {val}"
        items.append(key_val)

    return "\n".join(items)
