" functions for arranging foxed wifth columns"


def single_value_columns(values, print_width, col_width=40, is_truncated=False):
    """
    Display columns of single value.
    """

    output = ""
    line_nr = 0

    # Remove blank values
    values = [f"{value}" for value in values if value]

    # The number of columns
    col_count = print_width // col_width

    # The length of each column
    col_length = len(values) // col_count

    # Remove enough values so the total is dividable by column count
    values = values[: len(values) - (len(values) % col_count)]

    # Truncate the values that are longer than the column width.
    values = [f"{value[:col_width - 5]}..." if len(value) > col_width - 2 else value for value in values]

    # Pad each value so it fits the column width
    values = [f"{value:<{col_width}}" for value in values]

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


def key_val_columns(dict, print_width: int, col_width=40, ignore_keys=[], indent=""):
    """
    Display a dictionary's values in columns.
    """

    output = ""
    line_nr = 0
    gap = 4

    # The number of columns
    col_count = print_width // col_width

    # Remove blank values
    # for key, val in dict.items():

    # Gap between the columns.
    col_width = col_width - gap

    items = []

    for key, val in dict.items():
        key_len = len(str(key))
        val_len = len(str(val))

        if key_len + val_len + 1 > col_width:
            spacer_key = "<soft>" + ("." * (col_width - key_len - 1)) + "</soft>"
            spacer_val = "<soft>" + ("." * (col_width - val_len)) + "</soft>"
            val = f"{val[:col_width - 3]}..." if val_len > col_width else val
            key_val_1 = f"<cyan>{key}</cyan>:{spacer_key}"
            key_val_2 = f"{spacer_val}{val}"
            items.append(key_val_1 + (" ") * gap)
            items.append(key_val_2 + (" ") * gap)

        else:
            avail_key_width = col_width - val_len - 2
            if key_len > avail_key_width:
                key = key[: avail_key_width - 5] + "..."
            key_len = len(str(key))
            spacer = "<soft>" + ("." * (col_width - key_len - val_len - 1)) + "</soft>"
            key_val = f"<cyan>{key}:</cyan>{spacer}{val}"
            items.append(key_val + (" ") * gap)

    # The length of each column
    col_length = len(items) // col_count

    # Reorder the items to the print top to bottom, left to right
    # instead of left to right top to bottom
    cols = []

    for i in range(col_count):
        if i % col_count == 0:
            cols.append([])
        cols[cols[-1]].append(items[i])

    print(cols)

    # Compile the output
    for item in items:
        if line_nr % col_count == 0:
            output = output + "\n"
        output = output + item
        line_nr = line_nr + 1 % col_length

    return output


def process_list(values, cli_width, display_width: int, exclusions=[], indent=""):
    return_string = ""
    for item in values:
        if isinstance(item, dict):
            return_string = (
                return_string + "\n" + key_val_columns(item, cli_width, display_width, exclusions, "  " + indent)
            )
        else:
            return_string = return_string + "\n" + indent + "- " + str(item)
    return return_string
