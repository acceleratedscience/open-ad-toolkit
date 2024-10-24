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
    "This file is auto-generated.\n"
    "To update it, consult instructions:\n"
    "https://github.com/acceleratedscience/open-ad-toolkit/tree/main/docs\n\n"
    "-->"
)
DO_NOT_EDIT_PYPI = (
    "<!--\n\n"
    "DO NOT EDIT\n"
    "-----------\n"
    "This file is auto-generated with modified links for PyPI.\n"
    "To update it, consult instructions:\n"
    "https://github.com/acceleratedscience/open-ad-toolkit/tree/main/docs\n\n"
    "-->"
)

# endregion

############################################################
# region - README.md (GitHub)


# Update the README.md file with OpenAD description.
def update_github_readme_md(filename="README.md"):
    output_text(f"<h1>Updating <yellow>{filename}</yellow> with OpenAD description</h1>", pad_top=2)

    # Read README.md file content
    readme_md, err_msg = open_file(filename, return_err=True)
    if not readme_md:
        output_text(FLAG_ERROR, pad_top=1)
        output_error(err_msg)
        return

    # Read description source file content
    description_txt, err_msg = open_file("docs/source/description.txt", return_err=True)
    if not description_txt:
        output_text(FLAG_ERROR, pad_top=1)
        output_error(err_msg)
        return

    # Insert description
    readme_md_1 = readme_md.split("<!-- description -->")[0]
    readme_md_2 = readme_md.split("<!-- /description -->")[1]
    readme_md = readme_md_1 + "<!-- description -->\n" + description_txt + "\n<!-- /description -->" + readme_md_2

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
# region - README_plugins.md (GitHub)


# Update the README_plugins.md file with the about_plugin description
def update_github_readme_plugin_md(filename="README_plugin.md"):
    output_text(f"<h1>Updating <yellow>{filename}</yellow> with about_plugin description</h1>", pad_top=2)

    # Read README.md file content
    readme_plugin_md, err_msg = open_file(filename, return_err=True)
    if not readme_plugin_md:
        output_text(FLAG_ERROR, pad_top=1)
        output_error(err_msg)
        return

    # Read about_plugin source file content
    about_plugin_txt, err_msg = open_file("docs/source/about_plugin.txt", return_err=True)
    if not about_plugin_txt:
        output_text(FLAG_ERROR, pad_top=1)
        output_error(err_msg)
        return

    # Insert about_plugin text
    readme_md_1 = readme_plugin_md.split("<!-- about_plugin -->")[0]
    readme_md_2 = readme_plugin_md.split("<!-- /about_plugin -->")[1]
    readme_plugin_md = (
        readme_md_1 + "<!-- about_plugin -->\n" + about_plugin_txt + "\n<!-- /about_plugin -->" + readme_md_2
    )

    # Write to output file
    success, err_msg = write_file(filename, readme_plugin_md, return_err=True)
    if success:
        output_text(FLAG_SUCCESS)
        output_text(f"<soft>Updated</soft> <reset>{filename}</reset>")
    else:
        output_text(FLAG_ERROR)
        output_error(err_msg, pad=0)


# endregion

############################################################
# region - README-PYPI.md (PyPi)


# Create the README--PYPI.md for PyPI package page
# - - -
# A modified version of the README.md file with links pointing
# to the documentation website instead of other readme files.
def render_pypi_readme_md(filename="README--PYPI.md"):
    output_text(f"<h1>Generating <yellow>{filename}</yellow></h1>", pad_top=2)

    # Read README.md file content
    readme_md, err_msg = open_file("README.md", return_err=True)
    if not readme_md:
        output_text(FLAG_ERROR, pad_top=1)
        output_error(err_msg)
        return

    # Remove comments from README content
    readme_md = re.sub(r"<!--.*?-->\n?", "", readme_md, flags=re.DOTALL)

    # Adjust the links to play nice with just-the-docs
    readme_md = re.sub(
        r"^\[(.*?)\]: README_(.*?).md", lambda m: f"[{m.group(1)}]({m.group(2)}.html)", readme_md, flags=re.MULTILINE
    )

    # Insert DO NOT EDIT comment
    readme_md = DO_NOT_EDIT_PYPI + "\n\n" + readme_md

    # Write to output file
    success, err_msg = write_file(filename, readme_md, return_err=True)
    if success:
        output_text(FLAG_SUCCESS)
        output_text(f"<soft>Exported to</soft> <reset>/docs/output/markdown/{filename}</reset>")
    else:
        output_text(FLAG_ERROR)
        output_error(err_msg, pad=0)


# endregion

############################################################
# region - docs pages


# Convert all main pages
def render_docs_pages():
    for page in [
        "index.md",
        "installation.md",
        "getting-started.md",
        "models-service.md",
        "plugins.md",
        "ai-assistant.md",
        "developers.md",
    ]:
        _render_docs_page(page)


# Convert README_xxx.md page into a xxx.md page for the docs website
def _render_docs_page(filename):
    output_text(f"<h1>Generating <yellow>{filename}</yellow></h1>", pad_top=2)

    # Read xxx.md input file content
    input_md, err_msg = open_file(f"docs/input/{filename}", return_err=True)
    if not input_md:
        output_text(FLAG_ERROR, pad_top=1)
        output_error(err_msg)
        return

    # Read README_xxx.md base file content
    base_md_path = "README.md" if filename == "index.md" else f"README_{filename}"
    base_md, err_msg = open_file(base_md_path, return_err=True)
    if not base_md:
        output_text(FLAG_ERROR, pad_top=1)
        output_error(err_msg)
        return

    # Remove navigation links (index.md only)
    base_md = re.sub(r"<!-- navigation -->.*?<!-- /navigation -->\n*?<br>", "", base_md, flags=re.DOTALL)

    # Remove BACK links on top (all pages except index.md)
    base_md = re.sub(r"^<sub>.+?</sub>\n*?", "", base_md)

    # Remove comments
    base_md = re.sub(r"<!--.*?-->\n?", "", base_md, flags=re.DOTALL)

    # Remove superfloous linebreaks
    base_md = re.sub(r"\n\s*\n", "\n\n", base_md)

    # Trim space from start and end of file
    base_md = base_md.strip()

    # Adjust template links to play nice with just-the-docs
    base_md = re.sub(
        r"^\[(.*?)\]: README_(.*?).md",
        lambda m: f"[{m.group(1)}]: {m.group(2)}.html",
        base_md,
        flags=re.MULTILINE,
    )

    # Adjust regular links to play nice with just-the-docs
    base_md = re.sub(r"\[(.*?)\]\(README_(.*?).md\)", lambda m: f"[{m.group(1)}]({m.group(2)}.html)", base_md)

    # Insert DO NOT EDIT comment
    input_md = re.sub(r"{{DO_NOT_EDIT}}", DO_NOT_EDIT, input_md)

    # Insert base file content
    input_md = re.sub(r"{{CONTENT}}", base_md, input_md)

    # Write to output file
    success, err_msg = write_file(f"docs/output/markdown/{filename}", input_md, return_err=True)
    if success:
        output_text(FLAG_SUCCESS)
        output_text(f"<soft>Exported to</soft> <reset>/docs/output/markdown/{filename}</reset>")
    else:
        output_text(FLAG_ERROR)
        output_error(err_msg, pad=0)


# endregion

############################################################
# region - base-concepts.md


# Generate the base-concepts.md file for the documentation website.
def render_base_concepts_md(filename="base-concepts.md"):
    output_text(f"<h1>Generating <yellow>{filename}</yellow></h1>", pad_top=2)

    # Read base-concepts.md input file content
    base_concepts_md, err_msg = open_file("docs/input/base-concepts.md", return_err=True)
    if not base_concepts_md:
        output_text(FLAG_ERROR, pad_top=1)
        output_error(err_msg)
        return

    # Read about_workspace.txt source file content
    about_workspace, err_msg = open_file("docs/source/about_workspace.txt", return_err=True)
    if not about_workspace:
        output_text(FLAG_ERROR, pad_top=1)
        output_error(err_msg)
        return

    # Read about_mws.txt source file content (molecule working set)
    about_mws, err_msg = open_file("docs/source/about_mws.txt", return_err=True)
    if not about_mws:
        output_text(FLAG_ERROR, pad_top=1)
        output_error(err_msg)
        return

    # Read about_plugin.txt source file content
    about_plugin, err_msg = open_file("docs/source/about_plugin.txt", return_err=True)
    if not about_plugin:
        output_text(FLAG_ERROR, pad_top=1)
        output_error(err_msg)
        return

    # Read about_context.txt source file content
    about_context, err_msg = open_file("docs/source/about_context.txt", return_err=True)
    if not about_context:
        output_text(FLAG_ERROR, pad_top=1)
        output_error(err_msg)
        return

    # Read about_run.txt source file content
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
    # Update README files
    update_github_readme_md()
    update_github_readme_plugin_md()

    # Generate README for PyPI
    render_pypi_readme_md()

    # Turn README files into pages for the documentation website
    render_docs_pages()

    # Generate additional bespoke pages for documentation website
    render_base_concepts_md("base-concepts.md")
    render_commands_md("commands.md")

    # Render additional files
    render_commands_csv("commands.csv")
    render_description_txt("llm_description.txt")

    # Move all generated markdown files to the documentation repo
    docs = []
    for filename in os.listdir(f"{REPO_PATH}/docs/output/markdown"):
        docs.append(filename)
    copy_docs(docs)
