import os
import yaml
from openad.helpers.output import output_text
from openad.helpers.locale import localize
from openad.plugins.style_parser import wrap_text
import openad.core.help as core_help


def assemble_plugin_metadata(plugin_dir, cmd_help_dicts):
    """Assemble plugin metadata from from metadata and description files."""

    metadata = {}

    # Load metadata from file
    plugin_metadata = {}
    try:
        metadata_file = os.path.join(plugin_dir, "plugin_metadata.yaml")
        with open(metadata_file, "r", encoding="utf-8") as f:
            plugin_metadata = yaml.safe_load(f)
    except Exception:  # pylint: disable=broad-except
        pass

    # Load plugin description from file, or from metadata
    plugin_description = ""
    try:
        description_file = os.path.join(plugin_dir, "plugin_description.txt")
        with open(description_file, "r", encoding="utf-8") as f:
            plugin_description = f.read()
    except FileNotFoundError:
        plugin_description = plugin_metadata.get("description", "<error>No plugin description found</error>")
    except Exception as e:  # pylint: disable=broad-except
        plugin_description = plugin_metadata.get(
            "description", f"<error>An error occurred reading the plugin's description file</error>\n<soft>{e}<soft>"
        )

    # Organize commands
    commands_organized = core_help.organize_commands(cmd_help_dicts)
    commands_organized = commands_organized.get("_plugins", {})
    commands_organized = commands_organized.get(plugin_metadata.get("name", ""), {})

    # Assemble
    metadata["name"] = plugin_metadata.get("name", "Unnamed Plugin")
    metadata["namespace"] = plugin_metadata.get("namespace", None)
    metadata["description"] = plugin_description
    metadata["commands"] = commands_organized
    metadata["author"] = plugin_metadata.get("author", "Unknown")
    metadata["version"] = plugin_metadata.get("version", "0")
    metadata["website"] = plugin_metadata.get("website", None)
    metadata["github"] = plugin_metadata.get("github", None)

    return metadata


def display_plugin_overview(plugin_metadata) -> str:
    """Display plugin splash screen with description, metadata and commands."""

    output = []

    plugin_name = plugin_metadata.get("name")
    plugin_description = wrap_text(plugin_metadata.get("description"))
    plugin_commands = plugin_metadata.get("commands")
    plugin_author = plugin_metadata.get("author")
    plugin_version = plugin_metadata.get("version")
    plugin_website = plugin_metadata.get("website")
    plugin_github = plugin_metadata.get("github")

    # Title and description
    output.append(f"<yellow><reverse> PLUGIN </reverse></yellow><reverse> {plugin_name} </reverse>")
    output.append(f"<yellow>{'-' * 80}</yellow>")
    output.append(f"<soft>v{plugin_version} / Author: {plugin_author}</soft>\n")
    output.append(plugin_description)

    # Website + GitHub URLs
    if plugin_website:
        output.append(f"\n<soft>Website:</soft> <link>{plugin_website}</link>")
    if plugin_github:
        output.append(f"<soft>GitHub:</soft>  <link>{plugin_github}</link>")

    # Commands
    output.append("\n\n\n<h1>Available Commands</h1>")
    for category, commands in plugin_commands.items():
        if category == "_namespace":
            continue
        output.append(f"\n{category}")
        for cmd, description in commands:
            output.append(f"<cmd>{cmd}</cmd>")
    output = "\n".join(output)
    return output_text(output, edge=True, pad=2)


def reorder_commands_by_category_index(plugin_commands: dict) -> list:
    """
    Reorder the commands by their index, per category.

    Parameters
    ----------
    plugin_commands : dict
        The plugin commands, unordered, eg:
        [
            { category: "Category A", command: "Foo", index: 1 },
            { category: "Category B", command: "Baz", index: 2 },
            { category: "Category A", command: "Bar", index: 0 },
            { category: "Category B", command: "Foo", index: 1 },
            { category: "Category A", command: "Baz", index: 2 },
            { category: "Category B", command: "Bar", index: 0 },
        ]

    Returns
    -------
    list
        The reordered plugin commands, eg:
        [
            { category: "Category A", command: "Bar", index: 0 },
            { category: "Category A", command: "Foo", index: 1 },
            { category: "Category A", command: "Baz", index: 2 },
            { category: "Category B", command: "Bar", index: 0 },
            { category: "Category B", command: "Foo", index: 1 },
            { category: "Category B", command: "Baz", index: 2 },
    """

    # Organize commands by category
    plugin_commands_organized = {}
    for Cmd in plugin_commands:
        if Cmd.category not in plugin_commands_organized:
            plugin_commands_organized[Cmd.category] = []
        plugin_commands_organized[Cmd.category].append(Cmd)

    # Reorder the commands by their index, per category
    plugin_commands_output = []
    for cat_cmds in plugin_commands_organized.values():
        cat_cmds.sort(key=lambda x: x.index)
        for Cmd in cat_cmds:
            plugin_commands_output.append(Cmd)

    return plugin_commands_output


def all_cmds_note(plugin_name, plugin_namespace):
    """
    Return a localized note explaining how to display all available commands for a plugin.

    Displayed at the bottom of an individual plugin command's help output.

    Parameters
    ----------
    plugin_name : str
        The name of the plugin, eg. "My Plugin"
    plugin_namespace : str
        The plugin's namespace, eg. "mp"
    """
    # fmt: off
    text = {
        "en": f"To see all available <reset>{plugin_name}</reset> commands, run <cmd>{plugin_namespace}</cmd>",
        "fr": f"Pour voir toutes les commandes de <reset>{plugin_name}</reset> disponibles, exécutez <cmd>{plugin_namespace}</cmd>",
        "ja": f"利用可能なすべての <reset>{plugin_name}</reset> コマンドを表示するには、<cmd>{plugin_namespace}</cmd> を実行します。",
    }
    # fmt: on

    text = localize(text)
    text = f"\n\n<reset><reverse> i </reverse></reset> <soft>{text}</soft>"
    return text
