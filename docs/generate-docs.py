############################################################
# region - setup

"""
This script generates the commands.md and installation.md files
for the just-the-docs documentation.

- commands.md --> Generated from the command help
- installation.md --> Adapted from the main README.md

To generate:

    python3 docs/generate-docs.py

Output:
    
    docs/markdown/commands.md
    docs/markdown/installation.md

After being regenerated, copy the files over to the documentation repo.
"""

import os
import re
import glob
import json

from openad.app.main import RUNCMD as cmd_pointer
from openad.app.global_var_lib import _all_toolkits
from openad.toolkit.toolkit_main import load_toolkit
from openad.plugins.style_parser import tags_to_markdown
from openad.helpers.output import msg, output_error

# Get the path of this python file's parent folder
REPO_PATH = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

# endregion

############################################################
# region - commands.md


# Loop through all commands and export them to a markdown file
# that is ready to be included in the just-the-docs documentation.
def render_commands_md(filename):
    output = []  # Markdown
    toc = []  # Table of content

    # Just-the-docs markdown context
    jtd_identifier = (
        "---",
        "title: Commands",
        "layout: home",
        "nav_order: 4",
        "---",
    )
    output.append("\n".join(jtd_identifier) + "\n")

    # Intro comment
    comment = (
        "DO NOT EDIT",
        "-----------",
        "This file auto-generated.",
        "To update it, see openad/docs/generate-docs.py",
    )
    comment = "\n".join(comment)
    output.append(f"<!--\n\n{comment}\n\n-->" + "\n")

    # Parse main commands
    output.append(f"## OpenAD\n")
    toc.append(_toc_link("OpenAD"))
    cmds = cmd_pointer.current_help.help_current
    cmds_organized = _organize(cmds)
    _compile_section(output, toc, cmds_organized)

    # Parse tookit commands
    for toolkit_name in _all_toolkits:
        output.append(f"## {toolkit_name}\n\n")
        toc.append(_toc_link(toolkit_name))
        success, toolkit = load_toolkit(toolkit_name, from_repo=True)
        toolkit_cmds = toolkit.methods_help
        toolkit_cmds_organized = _organize(toolkit_cmds)
        _compile_section(output, toc, toolkit_cmds_organized)

    # Write output to file to this python file's parent folder
    toc = "### Table of Contents\n" + "\n".join(toc) + "\n"
    output = output[:2] + [toc] + output[2:]
    output = "\n".join(output)
    with open(f"{REPO_PATH}/docs/markdown/{filename}", "w") as f:
        f.write(output)


# Organize commands of a single section by category.
def _organize(cmds, toolkit_name=None):
    commands_organized = {}

    # Organize commands by category.
    for cmd in cmds:
        # Get command string.
        cmd_str = cmd["command"]
        cmd_description = cmd["description"]

        if "parent" in cmd and cmd["parent"]:
            cmd_str = "  -> " + cmd_str

        # Get category.
        category = cmd["category"] if "category" in cmd else "Uncategorized"

        # Organize by category.
        if category in commands_organized:
            commands_organized[category].append((cmd_str, cmd_description))
        else:
            commands_organized[category] = [(cmd_str, cmd_description)]

    return commands_organized


# Compile all commands of a single section.
def _compile_section(output, toc, cmds_organized):
    output.append('<details markdown="block">')
    output.append("<summary>See commands</summary>\n")
    for category in cmds_organized:
        output.append(f"### {category}\n")
        toc.append(_toc_link(category, 1))
        for cmd_str, cmd_description in cmds_organized[category]:
            output.append(f"`{cmd_str.strip()}`{{: .cmd }}\n{_parse_description(cmd_description)}<br>\n")
        output.append("<br>\n")
    output.append("</details>\n")


# Prepare the command description for proper rendering.
def _parse_description(description):
    description = tags_to_markdown(description)

    # Style notes as blockquotes, and ensure they're always
    # followed by an empty line, to avoid the next lines to
    # be treated as part of the blockquote.
    description = re.sub(
        r"(\*\*Note:\*\*.+?)(\n{1,})",
        lambda match: f"  > {match.group(1)}\n\n"
        if len(match.group(2)) == 1
        else f"  > {match.group(1)}{match.group(2)}",
        description,
        flags=re.MULTILINE,
    )

    # description = description.splitlines()
    # description = "\n".join([line.strip() for line in description])
    return description.strip()


# Convert a title to a markdown
# link for the table of contents.
# Foo Bar --> #foo-bar
def _toc_link(title, level=0):
    dash = "  " * level + "- "
    return f"{dash}[{title}](#{title.replace(' ', '-').lower()})"


# endregion

############################################################
# region - installation.md


# Adapt the README.md to be repurposed as
# instalation page for just-the-docs.
def render_installation_md(filename):
    readme = ""
    try:
        with open("README.md", "r") as f:
            readme = f.read()
    except BaseException as err:
        output_error(msg("err_readme", err, split=True), cmd_pointer, return_val=True)
        return

    # Remove all comments.
    readme = re.sub(r"<!--.*?-->", "", readme, flags=re.DOTALL)

    # Remove whitespaces from empty lines
    # and superfloous linebreaks.
    readme = re.sub(r"\n\s*\n", "\n\n", readme)

    # Trim space from start and end of file.
    readme = readme.strip()

    # Remove header and replace with just-the-docs
    # page identifier, plus add intro comment.
    splitter = "## Quick Install"
    jtd_identifier = (
        "---",
        "title: Installation",
        "layout: home",
        "nav_order: 2",
        "---",
    )
    comment = (
        "DO NOT EDIT",
        "-----------",
        "This file auto-generated from the main OpenAD README.md",
        "To update it, edit the main README.md and then regenerate this file.",
        "For instructions, see openad/docs/generate-docs.py",
    )
    comment = "\n".join(comment)
    comment = f"<!--\n\n{comment}\n\n-->" + "\n"
    readme = "\n".join(jtd_identifier) + "\n\n" + comment + "\n\n" + splitter + readme.split(splitter)[1]

    # Write to file
    with open(f"{REPO_PATH}/docs/markdown/{filename}", "w") as f:
        f.write(readme)


# endregion

############################################################

if __name__ == "__main__":
    render_commands_md("commands.md")
    render_installation_md("installation.md")
