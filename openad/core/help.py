""" The help Module"""

import re
import os
import webbrowser
from openad.helpers.files import open_file
from openad.helpers.general import singular, is_toolkit_installed, get_print_width, get_locale
from openad.helpers.output import output_error
from openad.helpers.output_msgs import msg
from openad.app.global_var_lib import GLOBAL_SETTINGS

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
    name: str,  # Name of the comand - used for ...?
    command: str,  # Command structure, used for help, docs, training
    description: str,  # Description of the command, used for help, docs, training
    note: str = None,  # Additional note to the command, only used in help (eg. To learn more about runs, run `run ?`)
    url: str = None,  # Currently not used - URL to the documentation of the command?
    category: str = "Uncategorized",  # Category used to organize the commands in help & docs
    parent: str = None,  # Parent command, only relevant for follow-up commands like `result open`
):
    """Create a help dictionary"""
    return {
        "category": category,
        "name": name,
        "command": command,
        "description": description,
        "note": note,
        "url": url,
        "parent": parent,
    }


def help_dict_create_v2(
    category: str,
    command: str,
    description_file: str = None,
    description: str = None,
    note: str | dict = None,
):
    """
    Create a help dictionary

    Parameters
    ----------
    category: str
        Category used to organize the commands in help & docs
    command: str
        Command structure, used for help, docs, training
    description: str, optional (not recommended)
        Option A: One-line description of the command, eg. "Say hello to the world"
    description_file: str, optional (recommended)
        option B: Command dir + filename of the description .txt file, eg. "hello_world/description.txt"
        If localized versions of the file exist, they will be prioritized based on the user's
        locale setting, eg. "hello_world/description_fr.txt".
    note: str | dict, optional
        Optional bottom note to the command description, eg. "To learn more about xyz, run `xyz ?`"
        Can be a dict with localized notes, eg. {"en": "Note in English", "fr": "Note en français"}
    """

    description = description or _try_localize_description_file(description_file)
    note = _try_localize_note_text(note)

    return {
        "category": category,
        "command": command,
        "description": description,
        "note": note,
    }


# Prioritize localized description files if available
def _try_localize_description_file(path):
    path_localized = _localize_path(path)
    if os.path.isfile(path_localized):
        path = path_localized

    description = open_file(path)
    return description


# Localize the description filename with the
# locale settings from your terminal.
def _localize_path(path: str):
    lang = get_locale("lang")
    if lang:
        path = path.replace(".txt", f"_{lang}.txt")
    return path


def _try_localize_note_text(note: str):
    if not note:
        return None

    # Simple string note
    elif isinstance(note, str):
        return note

    # Localized note
    elif isinstance(note, dict):
        lang = get_locale("lang")
        if lang and lang in note:
            note = note[lang]
        else:
            note = note.get("en", None)

    return note


def all_commands(
    available_commands: list,
    toolkit_name: str = None,
    toolkit_current: object = None,
    no_title: bool = False,
    cmd_pointer: object = None,
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
        output = []
        if not no_title:
            output.append(f'<h1>Available Commands - {toolkit_name if toolkit_name else "Main"}</h1>')
        if toolkit_name and not is_toolkit_installed(toolkit_name, cmd_pointer):
            err_msg = output_error(
                msg("fail_toolkit_not_installed", toolkit_name),
                return_val=True,
                jup_return_format="markdown_data",
                nowrap=True,
            )
            output.append(err_msg)
        elif len(commands_organized):
            if not no_title:
                output.append("")
            for category, available_commands in commands_organized.items():
                if no_title:
                    output.append(f"<h1>{category} Commands</h1>")
                else:
                    output.append(f"{category}:")
                for command_str in available_commands:
                    output.append(f"<cmd>{command_str}</cmd>")
                if not no_title:
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


def queried_commands(matching_commands: object, inp: str = None, starts_with_only: bool = False):
    """
    Return a styles list with all commands matching the query.

    Command: `<string> ?` or `? <string>`
    """
    inp_singular = singular(inp)
    if starts_with_only:
        output = [f'<yellow>Commands starting with "{inp}"</yellow>']
    else:
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
            command_str = re.sub(
                rf"(?<!<){re.escape(inp)}(s?)(?![^<>]*?>)", rf"<underline>{re.escape(inp)}\1</underline>", command_str
            )
        else:
            # String match --> only underline the matching string.

            command_str = re.sub(
                rf"(?<!<){re.escape(inp)}(?![^<>]*?>)", rf"<underline>{re.escape(inp)}</underline>", command_str
            )
        # Work around for escape character compensation in re
        command_str = command_str.replace("\ ", " ")
        command_str = command_str.replace("\[", "[")
        command_str = command_str.replace("\]", "]")
        command_str = command_str.replace("\<", "<")
        command_str = command_str.replace("\>", ">")

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
    print_width = get_print_width()

    # Style command.
    if GLOBAL_SETTINGS["display"] == "notebook":
        command_str = f"<cmd>{command_str}</cmd>"
        description = command["description"]
        note = f'\n<soft>{command["note"]}</soft>' if "note" in command and command["note"] is not None else None
    else:
        command_str = style(f"<cmd>{command_str}</cmd>", width=print_width)
        description = style(command["description"], width=print_width)
        note = (
            style(f'\n<soft>{command["note"]}</soft>', width=print_width)
            if "note" in command and command["note"] is not None
            else None
        )

    # Separator
    sep_len = min(a_len(command_str), print_width)
    sep = "<soft>" + sep_len * "-" + "</soft>"

    # Style description
    output = [command_str, sep, description]
    if note is not None:
        output.append(note)
    return "\n".join(output)


# Display advanced help
def advanced_help():
    """Call advanced Help"""

    # @future: build interactive command help using blessed
    webbrowser.open("https://acceleratedscience.github.io/openad-docs/commands.html")


class OpenadHelp:
    """OpenAD help Class"""

    help_orig = []
    help_current = []
    help_model_services = []
    help_plugins = []

    def add_help(self, help_dict):
        for i in help_dict:
            self.help_current.append(i)

    def reset_help(self):
        self.help_current = self.help_orig.copy()
        self.help_current.extend(self.help_model_services)
        self.help_current.extend(self.help_plugins)
