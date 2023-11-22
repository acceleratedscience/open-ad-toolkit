""" The help Module"""
import re
import shutil
from openad.helpers.general import singular, is_toolkit_installed
from openad.helpers.output import msg, output_text, output_error

# Importing our own plugins.
# This is temporary until every plugin is available as a public pypi package.
from openad.plugins.style_parser import style, a_len

# line_break = "<br>"
# heading = "<h3>"
# heading_end = "</h3>"
# bold_begin = "<b>"
# bold_end = "</b>"
# italic_begin = "<i>"
# italic_end = "</i>"
# spacing = " "
# six_spacing = "      "
# bullet = " - "


# Create the help dictionary object for a command.
def help_dict_create(
    name: str, command: str, description: str, url: str = None, category: str = "Uncategorized", parent: str = None
):
    """Create a help dictionary"""
    return {
        "category": category,
        "name": name,
        "command": command,
        "description": description,
        "url": url,
        "parent": parent,
    }


def all_commands(
    available_commands: list, toolkit_name: str = None, toolkit_current: object = None, cmd_pointer: object = None
):
    """
    Return xml string listing all available commands organized by category.

    Command: `?`
    """

    def _compile(available_commands, toolkit_name=None):
        commands_organized = {}

        # Organize commands by category.
        for command in available_commands:
            # Get command string.
            command_str = command["command"]
            if "parent" in command and command["parent"]:
                command_str = "  -> " + command_str

            # Get category.
            category = command["category"] if "category" in command else "Uncategorized"

            # Organize by category.
            if category in commands_organized:
                commands_organized[category].append(command_str)
            else:
                commands_organized[category] = [command_str]

        # Compile output.
        output = [f'<h1>Available Commands - {toolkit_name if toolkit_name else "Main"}</h1>']
        if toolkit_name and not is_toolkit_installed(toolkit_name, cmd_pointer):
            err_msg = output_error(
                msg("fail_toolkit_not_installed", toolkit_name, split=True), cmd_pointer, return_val=True, nowrap=True
            )
            output.append(err_msg)
        elif len(commands_organized):
            output.append("")
            for category, available_commands in commands_organized.items():
                output.append(f"{category}:")
                for command_str in available_commands:
                    output.append(f"<cmd>{command_str}</cmd>")
                output.append("")

            if toolkit_name:
                output.append(
                    f"<reverse> i </reverse> <soft>To learn more about the {toolkit_name} toolkit, run <cmd>{toolkit_name.lower()}</cmd>.</soft>"
                )
        else:
            output.append("<error>No commands found.</error>")

        return "\n".join(output)

    #
    #

    # Only list commands for a specific toolkit.
    if toolkit_name:
        return _compile(available_commands, toolkit_name)

    # List all commands, including selected toolkit.
    else:
        main_commands = _compile(available_commands)
        if toolkit_current:
            toolkit_commands = "\n\n\n" + _compile(toolkit_current.methods_help, toolkit_current.toolkit_name)
        else:
            toolkit_commands = ""
        return main_commands + toolkit_commands


def queried_commands(matching_commands: object, inp: str = None):
    """
    Return a styles list with all commands matching the query.

    Command: `<string> ?` or `? <string>`
    """

    inp_singular = singular(inp)
    output = [f'<yellow>Commands containing "{inp}"</yellow>']

    # First list commands that have an exact word match.
    if matching_commands["match_word"]:
        output = _append_matches(matching_commands["match_word"], inp_singular, output, match_word=True)

    # Then list commands that contain the string.
    if matching_commands["match_start"] or matching_commands["match_anywhere"]:
        # First list matches that start with the query.
        output = _append_matches(matching_commands["match_start"], inp_singular, output)

        # Then list matches that contain the query.
        output = _append_matches(matching_commands["match_anywhere"], inp_singular, output)

    return "\n".join(output)


def _append_matches(match_list, inp, output, match_word=False):
    for command in match_list:
        # command_str = command['command']
        # Display parent command in front of its follow up commands.
        command_str = (
            f'{command["parent"]} -> {command["command"]}'
            if "parent" in command and command["parent"]
            else command["command"]
        )

        # Underline matching string
        if match_word:
            # Exact word match --> underline both single and plural instances.
            command_str = re.sub(rf"(?<!<){inp}(s?)(?![^<>]*?>)", rf"<underline>{inp}\1</underline>", command_str)
        else:
            # String match --> only underline the matching string.
            command_str = re.sub(rf"(?<!<){inp}(?![^<>]*?>)", rf"<underline>{inp}</underline>", command_str)

        output.append(f"- <cmd>{command_str}</cmd>")
    return output


def command_details(command: list, cmd_pointer):
    """
    Return a single command with its description.

    Command: `<command> ?`
    """

    # Display parent command in front of its follow up commands.
    command_str = (
        f'{command["parent"]} -> {command["command"]}'
        if "parent" in command and command["parent"]
        else command["command"]
    )

    # Define optimal paragraph width.

    try:
        columns, lines = shutil.get_terminal_size()
    except Exception:
        columns = 155
    paragraph_width = min(columns - 5, 150)

    # Style command.
    if cmd_pointer.notebook_mode:
        command_str = f"<cmd>{command_str}</cmd>"
        description = command["description"]
    else:
        command_str = style(f"<cmd>{command_str}</cmd>", width=paragraph_width)
        description = style(command["description"], width=paragraph_width)

    # Separator
    sep_len = min(a_len(command_str), paragraph_width)
    sep = "<soft>" + sep_len * "-" + "</soft>"

    # Style description
    return "\n".join([command_str, sep, description])


# Display advanced help
def advanced_help():
    """Call advanced Help"""
    return "<warning>Advanced help is yet to be implemented.</warning>"


class OpenadHelp:
    """OpenAD help Class"""

    help_orig = []
    help_current = []

    def add_help(self, help_dict):
        for i in help_dict:
            self.help_current.append[i]

    def reset_help(self):
        self.help_current = self.help_orig.copy()
