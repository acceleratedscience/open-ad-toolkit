" functions for arranging foxed wifth columns"


def single_value_columns(values, print_width, col_width=40, is_truncated=False):
    """
    Display columns of single value.
    """

    output = ""
    line_nr = 0

    # The number of columns
    col_count = print_width // col_width

    # Remove blank values
    values = [f"{value}" for value in values if value]

    # The length of each column
    col_length = len(values) // col_count

    # Remove enough values so teh total is dividable by column count
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


def key_val_columns(values, cli_width, print_width: int, col_width=40, ignore_keys=[], indent=""):
    """
    Display a dictionary's values in columns.
    """

    output = ""
    i = 0

    # The number of columns
    col_count = print_width // col_width

    col_width = col_width - 4

    for key, val in values.items():
        if not key or not val:
            continue

        key_len = len(str(key))
        val_len = len(str(val))

        if key_len + val_len + 2 > col_width:
            key_val = f"<cyan>{key}</cyan>:\n{val:>{col_width}}"
        else:
            avail_key_width = col_width - val_len - 2
            if key_len > avail_key_width:
                key = key[: avail_key_width - 5] + "..."
            key_len = len(str(key))
            key = f"<cyan>{key}: </cyan>"
            key_val = f"{key}{val:%>{col_width - key_len - 2}}"

            # import re

            # Replace all consecutive % characters with a . then prepend the string with X and append with Z
            # key_val = re.sub(r"(%+)", r".", key_val)

        output = output + key_val + "\n"

    # for key, val in values.items():
    #     if key in ignore_keys:
    #         continue
    #     if val is None:
    #         continue
    #     buffer = print_width

    #     # if isinstance(value, dict):
    #     #     output = (
    #     #         output
    #     #         + "\n"
    #     #         + indent
    #     #         + f"<{key}:> "
    #     #         + "\n"
    #     #         + key_val_columns(value, cli_width, print_width, ignore_keys, indent="  " + indent)
    #     #     )
    #     #     continue

    #     # if isinstance(value, list):
    #     #     output = (
    #     #         output
    #     #         + "\n"
    #     #         + indent
    #     #         + f"<{key}:> "
    #     #         + process_list(value, cli_width, print_width, ignore_keys, indent="  " + indent)
    #     #     )
    #     #     continue

    #     value_string = f"<{key}:> {val}"

    #     while buffer < len(value_string) - 2:
    #         buffer = print_width + buffer
    #     if len(f"{value_string:<{print_width}}") + i < cli_width:
    #         output = output + indent + f"{value_string:<{buffer}}"
    #         i = i + len(f"{value_string:<{buffer}}")
    #     else:
    #         output = output + "\n" + indent + f"{value_string:<{buffer}}"
    #         i = len(f"{ value_string:<{buffer}}")

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
