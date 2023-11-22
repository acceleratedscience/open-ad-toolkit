import json

from openad.helpers.output import msg, output_text, output_error
from openad.helpers.ascii_type import ascii_type

# Importing our own plugins.
# This is temporary until every plugin is available as a public pypi package.
from openad.plugins.style_parser import style, wrap_text, strip_tags


def splash(toolkit_name=None, cmd_pointer=None, startup=False, raw=False):
    """Display the splash page for OpenAD or any of the toolkits."""

    toolkit_name = toolkit_name.upper() if toolkit_name else None

    # If no toolkit name is given, display the OpenAD welcome splash.
    main_splash = toolkit_name is None

    # Splash is minimized if toolkit is not active.
    toolkit_is_active = cmd_pointer and cmd_pointer.settings["context"] == toolkit_name

    # Toolkit splashes use full character.
    char = 0 if main_splash else 1

    # Load metadata from JSON file.
    import os

    if main_splash:
        json_file_path = f"/app/metadata.json"
    else:
        json_file_path = f"user_toolkits/{toolkit_name}/metadata.json"
    try:
        with open(os.path.dirname(os.path.abspath(__file__)) + "/../" + json_file_path) as json_file:
            data = json.load(json_file)
    except FileNotFoundError:
        return output_error(
            msg("fail_file_not_found", os.path.dirname(os.path.abspath(__file__)) + "/" + json_file_path, split=True)
        )

    # Make up for missing data.
    required_fields = "banner title intro author version commands".split()
    if main_splash:
        required_fields.remove("intro")
        required_fields.remove("author")
    for req_field in required_fields:
        if not data[req_field]:
            data[req_field] = "(not available)"

    # Compile output.
    output = ""

    # Startup message.
    # When you start the application with the context set to one of the toolkits,
    # you'll see the splash page for that toolkit. To make it more clear that this
    # is not the main splash page, we add a note on top clarifying this.
    if startup and not main_splash:
        output += f"<on_green> Your context is set to {toolkit_name}. </on_green>\n"
        output += "To see the main splash page, run <cmd>openad</cmd>.\n"
        output += f"To exit {toolkit_name}, run <cmd>unset context</cmd>.\n\n- - -\n\n"

    # Header
    if main_splash or toolkit_is_active:
        output += "<pre>" + ascii_type(data["banner"], reverse=not main_splash, char=char) + "</pre>" + "\n\n"
        output += f'<h1>{data["title"]}</h1>\n'
    else:
        output += f'<h1>{toolkit_name}: {data["title"]}</h1>\n'
    if main_splash:
        pass
        # output += '\n'
    else:
        output += f'<soft>v{data["version"]} / Author: {data["author"]}</soft>\n\n'

    # Intro
    if not main_splash:
        output += wrap_text(data["intro"]) + "\n\n"

    # Active
    if main_splash or toolkit_is_active:
        # Commands
        if not main_splash:
            data["commands"]["openad"] = "Display the main splash page."  # Add to commands for every toolkit.
        left_column_width, right_column_width = _calc_col_width(data["commands"])
        output += "\n" + _compile_commands(data["commands"], left_column_width, right_column_width)

    # Inactive
    elif cmd_pointer:
        # Installed
        if toolkit_name in cmd_pointer.settings["toolkits"]:
            output += msg("toolkit_installed", toolkit_name)

        # Not installed
        else:
            output += output_error(
                msg("fail_this_toolkit_not_installed", toolkit_name, split=True),
                cmd_pointer,
                return_val=True,
                jup_return_format="markdown_data",
                nowrap=True,
                pad=0,
            )

    if raw or cmd_pointer.notebook_mode:
        return output
    else:
        return style(output, pad=3, edge=True, nowrap=True)


# Calculate columns width.
def _calc_col_width(commands):
    left_column_width = 0
    gap = 3
    for key in commands:
        command_len = len(strip_tags(key)) + gap
        if command_len > left_column_width:
            left_column_width = command_len
    right_column_width = 80 - left_column_width
    return left_column_width, right_column_width


# Compile paragraph as array.
def _compile_commands(commands, left_column_width, right_column_width):
    output = []

    for key, val in commands.items():
        command = f"<cmd>{key}</cmd>"
        whitespace = (left_column_width - len(strip_tags(key))) * " "
        about = wrap_text(val, width=right_column_width)

        # Add necessary whitespace so the right column aligns nicely.
        def add_whitespace(i, line):
            if i == 0:
                return command + whitespace + line
            else:
                return left_column_width * " " + line

        lines = about.splitlines()
        about = list(map(lambda item: add_whitespace(*item), enumerate(lines)))

        # Merge multiple lines into one string.
        output.append("\n".join(about))

    return "\n\n".join(output)
