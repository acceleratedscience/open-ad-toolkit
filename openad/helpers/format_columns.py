" functions for arranging foxed wifth columns"


def single_value_columns(values, cli_width, display_width):
    """displays columns of single value"""
    return_string = ""
    i = 0
    for value in values:
        if len(str(value).strip()) == 0:
            continue

        if len(f"{value:<{display_width}}") + i < cli_width:
            return_string = return_string + f"{value:<{display_width}}"
            i = i + len(f"{value:<{display_width}}")
        else:
            return_string = return_string + f"\n{value:<{display_width}}"
            i = len(f"{value:<{display_width}}")
    return return_string


def name_and_value_columns(values, cli_width, display_width: int, exclusions=[], indent=""):
    """displays a nested dictionary"""
    return_string = ""
    i = 0

    for key, value in values.items():
        if key in exclusions:
            continue
        if value is None:
            continue
        buffer = display_width
        if isinstance(value, dict):
            return_string = (
                return_string
                + "\n"
                + indent
                + f"<{key}:> "
                + "\n"
                + name_and_value_columns(value, cli_width, display_width, exclusions, indent="  " + indent)
            )
            continue
        if isinstance(value, list):
            return_string = (
                return_string
                + "\n"
                + indent
                + f"<{key}:> "
                + process_list(value, cli_width, display_width, exclusions, indent="  " + indent)
            )
            continue
        value_string = "<{}:> {} ".format(key, str(value))

        while buffer < len(value_string) - 2:
            buffer = display_width + buffer
        if len(f"{value_string:<{display_width}}") + i < cli_width:
            return_string = return_string + indent + f"{value_string:<{buffer}}"
            i = i + len(f"{value_string:<{buffer}}")
        else:
            return_string = return_string + "\n" + indent + f"{value_string:<{buffer}}"
            i = len(f"{ value_string:<{buffer}}")

    return return_string


def process_list(values, cli_width, display_width: int, exclusions=[], indent=""):
    return_string = ""
    for item in values:
        if isinstance(item, dict):
            return_string = (
                return_string + "\n" + name_and_value_columns(item, cli_width, display_width, exclusions, "  " + indent)
            )
        else:
            return_string = return_string + "\n" + indent + "- " + str(item)
    return return_string
