""" The help Module"""

import re
import os
import webbrowser
from openad.helpers.general import singular, is_toolkit_installed, get_print_width
from openad.helpers.output import output_error
from openad.helpers.output_msgs import msg
from openad.helpers.locale import open_localized_file
from openad.helpers.plugins import all_cmds_note
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
        "commands": [command],
        "description": description,
        "note": note,
        "url": url,
        "parent": parent,
    }


def help_dict_create_v2(
    category: str,
    command: str | list,
    plugin_name: str = None,
    plugin_namespace: str = None,
    description_file: str = None,
    description: str = None,
):
    """
    Create a help dictionary.

    The help dictionaries are stored in cmd_pointer.current_help.help_plugins and are consumed by command_details()

    Parameters
    ----------
    category: str
        Category used to organize the commands in help & docs
    command: str | list
        Command string showing structure, used for help, docs, training
        Can be string of list of strings (to support multiple aliases)
    plugin_name: str, optional
        Name of the plugin this command belongs to, used to organize the commands in help & docs
    plugin_namespace: str, optional
        Namespace of the plugin this command belongs to, used to display help shortcut in help & docs
    description_file: str, optional
        option A: Command directory + filename of the description .txt file, eg. "hello_world/description.txt"
        If localized versions of the file exist, they will be prioritized based on the user's
        locale setting, eg. "hello_world/description_fr.txt".
    description: str, optional
        Option B: The actual desciption of the command, useful for one-line descriptions (eg. "Say hello to the world")
        or when the description needs to be parsed with variables (eg. f"Say hello to {name}"). When both description_file
        and description are provided, description will be prioritized.
    """

    description = description or open_localized_file(description_file)

    # Always store command string as a list, even if it's a single string.
    # Some commands have multiple aliases, so we need to support multiple strings.
    cmd_str_list = command if isinstance(command, list) else [command]

    return {
        "plugin_name": plugin_name,
        "plugin_namespace": plugin_namespace,
        "category": category,
        "commands": cmd_str_list,
        "description": description,
    }


def organize_commands(cmds):
    """
    Organize commands by category.

    When a command belongs to a plugin, it will be placed under the plugin's name
    under the '_plugins' key, with each plugin adhering to the same structure as
    the parent dictionary, except for the additinal '_namespace' key.

    Plugin commands are stored and organized separately, so they don't usually
    mix together with the main commands. When no plugin commands are present,
    the '_plugins' key will not be added.

    Input:
    [
        { "category": "Foo", "commands": ["foo a"], "description": <description> },
        { "category": "Foo", "commands": ["foo b"], "description": <description> },
        { "category": "Foo", "commands": ["foo c", "foo d"], "description": <description> },
        { "category": "Bar", "commands": ["bar a"], "description": <description> },
        { "category": "Bar", "commands": ["bar b"], "description": <description> },
        { "category": "Baz", "commands": ["baz a"], "description": <description>, "plugin_name": "Some Plugin", plugin_namespace: "sp" },
        { "category": "Baz", "commands": ["baz b"], "description": <description>, "plugin_name": "Some Plugin", plugin_namespace: "sp" },
    ]

    Output:
    {
        "Foo": [
            ("foo a", <description>),
            ("foo b", <description>),
            ("foo c", <description>),
            ("foo d", <description>),
        ],
        "Bar": [
            ("bar a", <description>),
            ("bar b", <description>),
        ],
        "_plugins": {
            "Some Plugin": {
                "_namespace": "sp",
                "Baz": [
                    ("baz a", <description>),
                    ("baz b", <description>),
                ],
            },
        },
    }
    """
    commands_organized = {}

    # Organize commands by category
    for cmd in cmds:
        # print('\n',cmd)
        # Get command description
        cmd_description = cmd.get("description")

        # Add parent command to child commands # Todo: remove
        if cmd.get("parent"):
            for cmd_str in cmd.get("commands"):
                cmd_str = "  -> " + cmd_str
        # Get category
        category = cmd.get("category") or "Uncategorized"
        # Get possible parent plugin
        plugin_name = cmd.get("plugin_name")

        # Organize plugin commands
        if plugin_name:
            if "_plugins" not in commands_organized:
                commands_organized["_plugins"] = {}
            if plugin_name not in commands_organized["_plugins"]:
                commands_organized["_plugins"][plugin_name] = {
                    "_namespace": cmd.get("plugin_namespace", plugin_name.lower())
                }
            for cmd_str in cmd.get("commands"):
                if category in commands_organized["_plugins"][plugin_name]:
                    commands_organized["_plugins"][plugin_name][category].append((cmd_str, cmd_description))
                else:
                    commands_organized["_plugins"][plugin_name][category] = [(cmd_str, cmd_description)]

        # Organize main commands
        else:
            # Organize by category
            if category in commands_organized:
                for cmd_str in cmd.get("commands"):
                    commands_organized[category].append((cmd_str, cmd_description))
            else:
                for cmd_str in cmd.get("commands"):
                    commands_organized[category] = [(cmd_str, cmd_description)]

    # Move plugins to the bottom
    if "_plugins" in commands_organized:
        plugins = commands_organized.pop("_plugins")
        commands_organized["_plugins"] = plugins

    return commands_organized


def all_commands(
    commands_organized: list,
    toolkit_name: str = None,
    plugin_name: str = None,
    toolkit_current: object = None,
    is_category: bool = False,
    cmd_pointer: object = None,
):
    """
    Return string for CLI printing, listing all commands that were passed, organized by category.

    This is used by the following commands:
        `?`
        `<category> ?` - eg. `molecule working set ?`
        `? <category>` - eg. `? molecule working set`
        `<toolkit> ?` - eg. `ds4sd ?`
        `? <toolkit>` - eg. `? ds4sd`
    """

    # Compile the help output.
    # - - -
    # Cycle through categories in the organized commands and add
    # commands lines and category titles to the output list.
    def _compile(commands_organized, toolkit_name=None, plugin_name=None):
        # Compile output.
        output = []

        # HEADER
        # - - - - - - - - - - - - - -

        # Category commands - no header
        if is_category:
            pass

        # Plugin commands - add header
        elif plugin_name:
            namespace = (
                commands_organized.pop("_namespace") if "_namespace" in commands_organized else plugin_name.lower()
            )
            output.append(f"<yellow><reverse> PLUGIN </reverse></yellow><reverse> {plugin_name} </reverse>\n")
            output.append(f"<soft>To learn more about this plugin, run <cmd>{namespace}</cmd></soft>\n")

        # Toolkit commands - add header
        elif toolkit_name:
            # Toolkit not installed
            if not is_toolkit_installed(toolkit_name, cmd_pointer):
                err_msg = output_error(
                    msg("fail_toolkit_not_installed", toolkit_name),
                    return_val=True,
                    jup_return_format="markdown_data",
                    nowrap=True,
                    pad_btm=1,
                )
                output.append(err_msg)

            #  Add toolkit header
            else:
                output.append(f"<h1>Available Commands - {toolkit_name}</h1>\n")

        # Main commands - add header
        else:
            output.append("<h1>Available Commands - Main</h1>\n")

        # COMMANDS
        # - - - - - - - - - - - - - -

        # Add commands
        if commands_organized:
            output = output + _add_cmds_by_category(commands_organized, is_plugin=bool(plugin_name))

            if toolkit_name:
                output.append(
                    f"\n<reverse> i </reverse> <soft>To learn more about the {toolkit_name} toolkit, run <cmd>{toolkit_name.lower()}</cmd>.</soft>"
                )
        else:
            output.append("<error>No commands found.</error>")

        return "\n".join(output)

    # Cycle through individual commands of a category and add them to the output list
    def _add_cmds_by_category(commands_organized: list, is_plugin: bool = False):
        output = []
        edge = "<soft>|</soft>    " if is_plugin and not is_category else ""

        for i, (category, available_commands) in enumerate(commands_organized.items()):
            # Ignore _plugins key, this is not a category but
            # a container for plugins - see organize_commands()
            if category == "_plugins":
                continue

            # First line blank, edge after that
            if i > 0:
                output.append(edge)

            # Category title
            if category != "Uncategorized":
                output.append(f"{edge}{category}")

            # Commands
            for cmd_str, _description in available_commands:
                for line in cmd_str.splitlines():
                    output.append(f"{edge}<cmd>{line}</cmd>")

            # Gap between the next category
            if not category:
                last_in_list = len(commands_organized.items()) - 1 == i
                output.append("" if last_in_list else edge)

        return output

    #
    #

    # Only list commands for a specific plugin.
    if plugin_name:
        return _compile(commands_organized, plugin_name=plugin_name)

    # Only list commands for a specific toolkit.
    if toolkit_name:
        return _compile(commands_organized, toolkit_name=toolkit_name)

    # List all commands, including selected toolkit.
    else:
        main_commands = _compile(commands_organized)
        toolkit_commands = ""
        plugin_commands = ""

        if not is_category:
            # List plugin commands
            all_plugin_commands_organized = organize_commands(cmd_pointer.current_help.help_plugins).get("_plugins", {})
            for plugin_name, plugin_commands_organized in all_plugin_commands_organized.items():
                plugin_commands = (
                    plugin_commands + "\n\n\n\n" + _compile(plugin_commands_organized, plugin_name=plugin_name)
                )

            # List toolkit commands
            if toolkit_current:
                toolkit_commands_organized = organize_commands(toolkit_current.methods_help)
                toolkit_commands = "\n\n\n\n" + _compile(
                    toolkit_commands_organized, toolkit_name=toolkit_current.toolkit_name
                )
            else:
                toolkit_commands = ""

        return main_commands + toolkit_commands + plugin_commands


def queried_commands(matching_commands: object, inp: str = None, starts_with_only: bool = False):
    """
    Return a styled list with all commands matching the query.

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
    for cmd in match_list:
        # Style command(s)
        for cmd_str in cmd.get("commands"):
            # Display parent command in front of its follow up commands.
            cmd_str = f'{cmd["parent"]} -> {cmd_str}' if cmd.get("parent") else cmd_str

            # Underline matching string
            if match_word:
                # Exact word match --> underline both single and plural instances.
                cmd_str = re.sub(
                    rf"(?<!<){re.escape(inp)}(s?)(?![^<>]*?>)", rf"<underline>{re.escape(inp)}\1</underline>", cmd_str
                )
            else:
                # String match --> only underline the matching string.
                cmd_str = re.sub(
                    rf"(?<!<){re.escape(inp)}(?![^<>]*?>)", rf"<underline>{re.escape(inp)}</underline>", cmd_str
                )

            # Work around for escape character compensation in re
            cmd_str = cmd_str.replace("\ ", " ")
            cmd_str = cmd_str.replace("\[", "[")
            cmd_str = cmd_str.replace("\]", "]")
            cmd_str = cmd_str.replace("\<", "<")
            cmd_str = cmd_str.replace("\>", ">")

            output.append(f"- <cmd>{cmd_str}</cmd>")
    return output


def command_details(cmd: list):
    """
    Return a single command with its description.

    Command: `<command> ?`
    """

    # Define optimal paragraph width
    print_width = get_print_width()

    # Style command(s)
    cmd_str_list_styled = []
    sep_len = 0
    for cmd_str in cmd.get("commands"):
        # Display parent command in front of its follow up commands
        cmd_str = f'{cmd["parent"]} -> {cmd_str}' if cmd.get("parent") else cmd_str
        sep_len = min(max(a_len(cmd_str), sep_len), print_width)

        # Style command
        if GLOBAL_SETTINGS["display"] == "notebook":
            cmd_str_list_styled.append(f"<cmd>{cmd_str}</cmd>")
        else:
            cmd_str_list_styled.append(style(f"<cmd>{cmd_str}</cmd>", width=print_width))

    # Bottom note
    note = None
    if "plugin_name" in cmd and "plugin_namespace" in cmd:
        note = all_cmds_note(cmd.get("plugin_name"), cmd.get("plugin_namespace"))

    # Style description and note
    if GLOBAL_SETTINGS["display"] == "notebook":
        description = cmd["description"]

    else:
        description = style(cmd["description"], width=print_width)
        note = style(note, width=print_width) if note else None

    # Separator
    sep = "<soft>" + sep_len * "-" + "</soft>"

    # Style description
    output = ["\n".join(cmd_str_list_styled), sep, description]
    if note:
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
        # self.help_current.extend(self.help_plugins)
