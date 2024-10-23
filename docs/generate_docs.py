"""
Generate documentation files for the OpenAD Toolkit.
For more information, consult the README.md file in the docs folder.

python3 docs/generate_docs.py

"""

############################################################
# region - setup

import os
import re
import sys
import pyperclip

# Add the root directory to the sys.path
root_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
if str(root_dir) not in sys.path:
    sys.path.append(root_dir)
# for path in sys.path:
#     print("*", path)

from copy_docs import copy_docs  # This resolves when running the script directly
from openad.app.main import RUNCMD as cmd_pointer
from openad.app.global_var_lib import _all_toolkits
from openad.toolkit.toolkit_main import load_toolkit
from openad.plugins.style_parser import tags_to_markdown
from openad.helpers.output import output_error, output_text, output_success
from openad.helpers.output_msgs import msg
from openad.helpers.files import open_file, write_file

# Get the repo path, this python file's parent folder.
REPO_PATH = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
FLAG_SUCCESS = f"<on_green> SUCCESS </on_green>"
FLAG_ERROR = f"<on_red> FAILED </on_red>"
DO_NOT_EDIT = (
    "<!--\n\n"
    "DO NOT EDIT\n"
    "-----------\n"
    "This file auto-generated.\n"
    "To update it, consult instructions:\n"
    "https://github.com/acceleratedscience/open-ad-toolkit/tree/main/docs\n\n"
    "-->"
)

# endregion

############################################################
# region - README.md


# Replace the description in the main README.md file.
def update_readme_md(filename):
    output_text("<h1>Updating <yellow>README.md</yellow> with OpenAD description</h1>", pad_top=2)

    # Read README.md input content
    readme_md, err_msg = open_file("README.md", return_err=True)
    if not readme_md:
        output_text(FLAG_ERROR, pad_top=1)
        output_error(err_msg)
        return

    # Read description file input content
    description_txt, err_msg = open_file("docs/source/description.txt", return_err=True)
    if not description_txt:
        output_text(FLAG_ERROR, pad_top=1)
        output_error(err_msg)
        return

    # Insert description
    readme_md_1 = readme_md.split("<!-- description -->")[0]
    readme_md_2 = readme_md.split("<!-- /description -->")[1]
    readme_md = readme_md_1 + "<!-- description -->\n\n" + description_txt + "\n\n<!-- /description -->" + readme_md_2

    # Write to output file
    success, err_msg = write_file(filename, readme_md, return_err=True)
    if success:
        output_text(FLAG_SUCCESS)
        output_text(f"<soft>Updated</soft> <reset>{filename}</reset>")
    else:
        output_text(FLAG_ERROR)
        output_error(err_msg, pad=0)


# endregion

############################################################
# region - index.md


# Loop through all commands and export them to a markdown file
# that is ready to be included in the just-the-docs documentation.
def render_index_md(filename):
    output_text("<h1>Generating <yellow>index.md</yellow></h1>", pad_top=2)

    # Read index.md input content
    index_md, err_msg = open_file("docs/input/index.md", return_err=True)
    if not index_md:
        output_text(FLAG_ERROR, pad_top=1)
        output_error(err_msg)
        return

    # Read description file input content
    description_txt, err_msg = open_file("docs/source/description.txt", return_err=True)
    if not description_txt:
        output_text(FLAG_ERROR, pad_top=1)
        output_error(err_msg)
        return

    # Insert DO NOT EDIT comment
    index_md = re.sub(r"{{DO_NOT_EDIT}}", DO_NOT_EDIT, index_md, flags=re.DOTALL)

    # Insert description
    index_md = re.sub(r"{{DESCRIPTION}}", description_txt, index_md, flags=re.DOTALL)

    # Write to output file
    success, err_msg = write_file(f"docs/output/markdown/{filename}", index_md, return_err=True)
    if success:
        output_text(FLAG_SUCCESS)
        output_text(f"<soft>Exported to</soft> <reset>/docs/output/markdown/{filename}</reset>")
    else:
        output_text(FLAG_ERROR)
        output_error(err_msg, pad=0)


# endregion

############################################################
# region - installation.md


# Adapt the README.md to be repurposed as
# instalation page for just-the-docs.
def render_installation_md(filename):
    output_text("<h1>Generating <yellow>installation.md</yellow> based off of README.md</h1>", pad_top=2)

    # Read installation.md input content
    installation_md, err_msg = open_file("docs/input/installation.md", return_err=True)
    if not installation_md:
        output_text(FLAG_ERROR, pad_top=1)
        output_error(err_msg)
        return

    # Read README.md content
    readme, err_msg = open_file("README.md", return_err=True)
    if not readme:
        output_text(FLAG_ERROR, pad_top=1)
        output_error(err_msg)
        return

    # Remove all comments.
    readme = re.sub(r"<!--.*?-->", "", readme, flags=re.DOTALL)

    # Remove whitespaces from empty lines
    # and superfloous linebreaks.
    readme = re.sub(r"\n\s*\n", "\n\n", readme)

    # Trim space from start and end of file.
    readme = readme.strip()

    # Remove header
    splitter = "## Quick Install"
    readme = splitter + readme.split(splitter)[1]

    # Insert DO NOT EDIT comment
    installation_md = re.sub(r"{{DO_NOT_EDIT}}", DO_NOT_EDIT, installation_md, flags=re.DOTALL)

    # Insert description
    installation_md = re.sub(r"{{INSTALLATION}}", readme, installation_md, flags=re.DOTALL)

    # Write to file
    success, err_msg = write_file(f"{REPO_PATH}/docs/output/markdown/{filename}", installation_md, return_err=True)
    if success:
        output_text(FLAG_SUCCESS)
        output_text(f"<soft>Exported to</soft> <reset>/docs/output/markdown/{filename}</reset>")
    else:
        output_text(FLAG_ERROR)
        output_error(err_msg, pad=0)


# endregion

############################################################
# region - base-concepts.md


# Adapt the README.md to be repurposed as
# instalation page for just-the-docs.
def render_base_concepts_md(filename):
    output_text("<h1>Generating <yellow>base-concepts.md</yellow></h1>", pad_top=2)

    # Read base-concepts.md input content
    base_concepts_md, err_msg = open_file("docs/input/base-concepts.md", return_err=True)
    if not base_concepts_md:
        output_text(FLAG_ERROR, pad_top=1)
        output_error(err_msg)
        return

    # Read about_workspace.txt content
    about_workspace, err_msg = open_file("docs/source/about_workspace.txt", return_err=True)
    if not about_workspace:
        output_text(FLAG_ERROR, pad_top=1)
        output_error(err_msg)
        return

    # Read about_mws.txt content
    about_mws, err_msg = open_file("docs/source/about_mws.txt", return_err=True)
    if not about_mws:
        output_text(FLAG_ERROR, pad_top=1)
        output_error(err_msg)
        return

    # Read about_plugin.txt content
    about_plugin, err_msg = open_file("docs/source/about_plugin.txt", return_err=True)
    if not about_plugin:
        output_text(FLAG_ERROR, pad_top=1)
        output_error(err_msg)
        return

    # Read about_context.txt content
    about_context, err_msg = open_file("docs/source/about_context.txt", return_err=True)
    if not about_context:
        output_text(FLAG_ERROR, pad_top=1)
        output_error(err_msg)
        return

    # Read about_run.txt content
    about_run, err_msg = open_file("docs/source/about_run.txt", return_err=True)
    if not about_run:
        output_text(FLAG_ERROR, pad_top=1)
        output_error(err_msg)
        return

    # Insert DO NOT EDIT comment
    base_concepts_md = re.sub(r"{{DO_NOT_EDIT}}", DO_NOT_EDIT, base_concepts_md, flags=re.DOTALL)

    # Insert descriptions
    base_concepts_md = re.sub(r"{{ABOUT_WORKSPACE}}", about_workspace, base_concepts_md, flags=re.DOTALL)
    base_concepts_md = re.sub(r"{{ABOUT_MWS}}", about_mws, base_concepts_md, flags=re.DOTALL)
    base_concepts_md = re.sub(r"{{ABOUT_PLUGIN}}", about_plugin, base_concepts_md, flags=re.DOTALL)
    base_concepts_md = re.sub(r"{{ABOUT_CONTEXT}}", about_context, base_concepts_md, flags=re.DOTALL)
    base_concepts_md = re.sub(r"{{ABOUT_RUN}}", about_run, base_concepts_md, flags=re.DOTALL)

    # Write to file
    success, err_msg = write_file(f"{REPO_PATH}/docs/output/markdown/{filename}", base_concepts_md, return_err=True)
    if success:
        output_text(FLAG_SUCCESS)
        output_text(f"<soft>Exported to</soft> <reset>/docs/output/markdown/{filename}</reset>")
    else:
        output_text(FLAG_ERROR)
        output_error(err_msg, pad=0)


# endregion

############################################################
# region - commands.md


# Loop through all commands and export them to a markdown file
# that is ready to be included in the just-the-docs documentation.
def render_commands_md(filename):
    output_text("<h1>Generating <yellow>commands.md</yellow> from help</h1>", pad_top=2)

    toc = []  # Table of content
    commands = []  # Markdown

    # Parse main commands
    commands.append(f"## OpenAD\n")
    toc.append(_toc_link("OpenAD"))
    cmds = cmd_pointer.current_help.help_current
    cmds_organized = _organize(cmds)
    _compile_section(commands, toc, cmds_organized)

    # Parse tookit commands
    for toolkit_name in _all_toolkits:
        commands.append(f"## {toolkit_name}\n\n")
        toc.append(_toc_link(toolkit_name))
        success, toolkit = load_toolkit(toolkit_name, from_repo=True)
        if success:
            toolkit_cmds = toolkit.methods_help
            toolkit_cmds_organized = _organize(toolkit_cmds)
            _compile_section(commands, toc, toolkit_cmds_organized)

    # Compile table of contents
    toc = "### Table of Contents\n" + "\n".join(toc) + "\n"

    # Compile commands
    commands = "\n".join(commands)

    # Read commands.md input content
    commands_md, err_msg = open_file("docs/input/commands.md", return_err=True)
    if not commands_md:
        output_text(FLAG_ERROR, pad_top=1)
        output_error(err_msg)
        return

    # Insert DO NOT EDIT comment
    commands_md = re.sub(r"{{DO_NOT_EDIT}}", DO_NOT_EDIT, commands_md, flags=re.DOTALL)

    # Insert table of contents
    commands_md = re.sub(r"{{TOC}}", toc, commands_md, flags=re.DOTALL)

    # Insert commands
    commands_md = re.sub(r"{{COMMANDS}}", commands, commands_md, flags=re.DOTALL)

    # Write to file
    success, err_msg = write_file(f"{REPO_PATH}/docs/output/markdown/{filename}", commands_md, return_err=True)
    if success:
        output_text(FLAG_SUCCESS)
        output_text(f"<soft>Exported to</soft> <reset>/docs/output/markdown/{filename}</reset>")
    else:
        output_text(FLAG_ERROR)
        output_error(err_msg, pad=0)


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
        lambda match: (
            f"  > {match.group(1)}\n\n" if len(match.group(2)) == 1 else f"  > {match.group(1)}{match.group(2)}"
        ),
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
# region - commands.csv


# Loop through all commands and export them to a CSV file.
# This is not used for anything in particular, other than
# to have a list of all commands in a file which can be annotated.
def render_commands_csv(filename, delimiter=";"):
    output_text("<h1>Generating <yellow>commands.csv</yellow> from help</h1>", pad_top=2)
    output = [["Command", "Category"]]

    # Parse main commands
    cmds_main = cmd_pointer.current_help.help_current
    cmds_organized = _organize(cmds_main)

    # Parse tookit commands
    for toolkit_name in _all_toolkits:
        success, toolkit = load_toolkit(toolkit_name, from_repo=True)
        if success:
            toolkit_cmds = toolkit.methods_help
            toolkit_cmds_organized = _organize(toolkit_cmds)
            cmds_organized.update(toolkit_cmds_organized)

    # Add a row per command.
    for category, cmds in cmds_organized.items():
        for cmd in cmds:
            output.append([cmd[0], category])

    # Convert to CSV string
    output_str = "\n".join([f"{delimiter}".join(row) for row in output])

    # Convert to clipboard CSV string
    output_clipboard = "\n".join([f"\t".join(row) for row in output])
    pyperclip.copy(output_clipboard)

    # Write to file
    success, err_msg = write_file(f"{REPO_PATH}/docs/output/csv/{filename}", output_str, return_err=True)
    if success:
        output_text(FLAG_SUCCESS)
        output_text(f"<soft>Exported to</soft> <reset>/docs/output/csv/{filename}</reset>")
    else:
        output_text(FLAG_ERROR)
        output_error(err_msg, pad=0)
    output_success(msg("csv_to_clipboard"), pad=0)


# endregion

############################################################
# region - toolkits -> llm_description.txt


# Update commands in the llm_description.txt file per toolkit.
# Used as training data by the LLM for the "tell me" command.
# - - -
# Note: llm_description.txt needs to be set up with the toolkit
# LLM briefing set up and the following line will define where
# the commands are to be inserted - any text after this line
# will be overwritten:
# "The following commands are available for this toolkit:"
def render_description_txt(filename):
    output_text("<h1>Updating commands in <yellow>llm_description.txt</yellow> for all toolkits</h1>", pad_top=4)

    # Loop through all toolkits
    for toolkit_name in _all_toolkits:
        flag_toolkit = f"<on_white><black> {toolkit_name} </black></on_white>"
        # Load toolkit
        success, toolkit = load_toolkit(toolkit_name, from_repo=True)
        if not success:
            err_msg = toolkit
            output_text(flag_toolkit + FLAG_ERROR)
            output_error(msg("err_load_toolkit", toolkit_name), pad=0)
            continue

        toolkit_cmds = toolkit.methods_help
        toolkit_cmds_organized = _organize(toolkit_cmds)
        output = _compile_commands(toolkit_cmds_organized)

        # Load llm_description.txt
        file_path = f"{REPO_PATH}/openad/user_toolkits/{toolkit_name}/{filename}"
        description_txt, err_msg = open_file(file_path, return_err=True)
        if not description_txt:
            output_text(flag_toolkit + FLAG_ERROR)
            # output_error(msg("err_load_toolkit_description", toolkit_name), pad=0) # Maybe overkill
            output_error(err_msg, pad_btm=1)
            continue

        # Insert commands into llm_description.txt
        splitter = "The following commands are available for this toolkit:"
        if splitter not in description_txt:
            output_text(flag_toolkit + FLAG_ERROR)
            output_error(msg("err_invalid_description_txt", toolkit_name, splitter), pad_btm=1)
            continue
        description_txt = description_txt.split(splitter)[0] + splitter + "\n\n"
        description_txt += "\n".join(output)
        description_txt = description_txt.strip()

        # print(("----" * 50) + "\n" + description_txt + "\n" + ("----" * 50))

        # Write to file
        success, err_msg = write_file(file_path, description_txt, return_err=True)
        if success:
            output_text(flag_toolkit + FLAG_SUCCESS)
            output_text(
                f"<soft>Updated in</soft> <reset>/docs/openad/user_toolkits/<toolkit_name>/{filename}</reset>",
                pad_btm=1,
            )
        else:
            output_text(flag_toolkit + FLAG_ERROR)
            output_error(err_msg, pad_btm=1)

    output_text("", pad_btm=2)


# Compile all commands for a single toolkit's llm_description.txt.
def _compile_commands(cmds_organized):
    output = []
    for category in cmds_organized:
        output.append(category + ":")
        for cmd_str, cmd_description in cmds_organized[category]:
            # Add command
            output.append(f"\t`{cmd_str.strip()}`")

            # Add command description
            cmd_description = tags_to_markdown(cmd_description).strip()
            cmd_description = cmd_description.replace("<br>", "")
            cmd_description = cmd_description.splitlines()
            cmd_description = "\n\t\t".join([line.strip() for line in cmd_description])
        output.append("")

    return output


# endregion

############################################################

if __name__ == "__main__":
    # Update main README.md
    update_readme_md("README.md")

    # Render markdown files for documentation website
    render_index_md("index.md")
    render_installation_md("installation.md")
    render_base_concepts_md("base-concepts.md")
    render_commands_md("commands.md")
    render_commands_csv("commands.csv")
    render_description_txt("llm_description.txt")

    # Copy markdown files into the documentation repo
    copy_docs(["index.md", "installation.md", "base-concepts.md", "commands.md"])
