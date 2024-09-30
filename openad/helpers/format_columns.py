" functions for arranging foxed wifth columns"

import math

# from openad.app.global_var_lib import PRINT_WIDTH
PRINT_WIDTH = 150


def single_value_columns(values: list, print_width=PRINT_WIDTH, col_width=40, is_truncated=False, indent=0):
    """
    Display columns of single value per line.
    """
    print(11)

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
        if line_nr % col_count == 0:
            output = output + "\n"
        output = output + val
        line_nr = line_nr + 1 % col_length

    # Add elipsis to the end of each column if the output is truncated
    if is_truncated:
        output = output + "\n" + (f"{'<soft>...</soft>':<{col_width + 13}}" * col_count)

    return output


def key_val_columns(items_dict: dict, print_width: int = PRINT_WIDTH, col_width=40, ignore_keys=[], indent=0):
    """
    Display a dictionary's values in columns.
    """
    print(12)

    print_width = print_width - indent
    output = ""
    line_nr = 0
    gap = 4

    # The number of columns
    col_count = print_width // col_width

    # Gap between the columns.
    col_width = col_width - gap

    items = []
    for key, val in items_dict.items():
        if key in ignore_keys:
            continue

        key_len = len(str(key))
        val_len = len(str(val))

        # Replace blank values with '-'
        if val != 0 and not val:
            key = "<soft>" + key + "</soft>"
            val = "<soft>-</soft>"
            val_len = 1

        if key_len + 2 + val_len > col_width:
            spacer_key = "<soft>" + ("." * (col_width - key_len - 1)) + "</soft>"
            spacer_val = "<soft>" + ("." * (col_width - val_len)) + "</soft>"
            val = f"{val[:col_width - 3]}..." if val_len > col_width else val
            key_color = "soft" if val != 0 and not val else "cyan"
            key_val_1 = f"{(indent * ' ')}<{key_color}>{key}</{key_color}>:{spacer_key}"
            key_val_2 = f"{(indent * ' ')}{spacer_val}{val}"
            items.append(key_val_1 + (" ") * gap)
            items.append(key_val_2 + (" ") * gap)

        else:
            spacer = "<soft>" + ("." * (col_width - key_len - val_len - 1)) + "</soft>"
            key_val = f"{(indent * ' ')}<cyan>{key}:</cyan>{spacer}{val}"
            items.append(key_val + (" ") * gap)

    # The length of each column
    col_length = math.ceil(len(items) / col_count)

    # Reorder the items to the print top to bottom, left to right
    # instead of left to right top to bottom
    cols = []
    for j, item in enumerate(items):
        if j % col_length == 0:
            cols.append([])
        cols[-1].append(item)

    items_rearranged = []
    i = 0
    j = 0
    while len(items_rearranged) < len(items) and j < 50:
        for col in cols:
            # print("%\n", len(col), col, "\n\n")
            if j < len(col):
                items_rearranged.append(col[j])
        i = i + 1
        j = (j + 1) % len(col)

    items = items_rearranged

    # Compile the output
    for item in items:
        if line_nr % col_count == 0:
            output = output + "\n"
        output = output + item
        line_nr = line_nr + 1 % col_length

    return output


def pretty_key_val(items_dict: dict, print_width: int = PRINT_WIDTH, ignore_keys=[], indent=0):
    items = []
    for key, val in items_dict.items():
        val_len = len(str(val))
        val = f"{val[:print_width - len(key) - 2 - 3]}..." if val_len > print_width else val
        key_val = f"{(indent * ' ')}<cyan>{key}:</cyan> {val}"
        items.append(key_val)

    return "\n" + "\n".join(items)


# # Unused ?
# def process_list(values, cli_width, display_width: int, ignore_keys=[], indent=0):
#     return_string = ""
#     for item in values:
#         if isinstance(item, dict):
#             return_string = (
#                 return_string + "\n" + key_val_columns(item, cli_width, display_width, ignore_keys, "  " + indent * " ")
#             )
#         else:
#             return_string = return_string + "\n" + indent * " " + "- " + str(item)
#     return return_string
