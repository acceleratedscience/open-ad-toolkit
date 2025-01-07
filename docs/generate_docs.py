"""
Generate documentation files for the OpenAD Toolkit.
For more information, consult the README.md file in the docs folder.

python3 docs/generate_docs.py

"""

############################################################
# region - setup

import os
import re
import pyperclip

# Add the root directory to the sys.path
# root_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
# if str(root_dir) not in sys.path:
#     sys.path.append(root_dir)
# for path in sys.path:
#     print("*", path)

from copy_docs import copy_docs  # This resolves when running the script directly
from openad.app.main import RUNCMD as cmd_pointer
from openad.app.global_var_lib import _all_toolkits
from openad.core.help import organize_commands
from openad.toolkit.toolkit_main import load_toolkit
from openad.plugins.style_parser import tags_to_markdown
from openad.helpers.output import output_error, output_text, output_success
from openad.helpers.output_msgs import msg
from openad.helpers.files import open_file, write_file

# Get the repo path, this python file's parent folder.
# REPO_PATH = os.path.dirname(os.path.dirname(os.path.realpath(__file__))) %%
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
    description_txt, err_msg = open_file("openad/docs_src/description.txt", return_err=True)
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
        output_text(f"<soft>Updated</soft> <reset>/{filename}</reset>")
    else:
        output_text(FLAG_ERROR)
        output_error(err_msg, pad=0)


# endregion

############################################################
# region - README/plugins.md (GitHub)


# Update the README_plugins.md file with the about_plugin description
def update_github_readme_plugins_md(filename="plugins.md"):
    output_text(f"<h1>Updating <yellow>{filename}</yellow> with about_plugin description</h1>", pad_top=2)

    # Read README.md file content
    readme_plugin_md, err_msg = open_file(f"README/{filename}", return_err=True)
    if not readme_plugin_md:
        output_text(FLAG_ERROR, pad_top=1)
        output_error(err_msg)
        return

    # Read about_plugin source file content
    about_plugin_txt, err_msg = open_file("openad/docs_src/about_plugin.txt", return_err=True)
    if not about_plugin_txt:
        output_text(FLAG_ERROR, pad_top=1)
        output_error(err_msg)
        return

    # Add some extra formatting to the about_plugin text
    about_plugin_txt = re.sub(r"^Note: ", "> **Note:** ", about_plugin_txt, flags=re.MULTILINE)

    # Insert about_plugin text
    readme_md_1 = readme_plugin_md.split("<!-- about_plugin -->")[0]
    readme_md_2 = readme_plugin_md.split("<!-- /about_plugin -->")[1]
    readme_plugin_md = (
        readme_md_1 + "<!-- about_plugin -->\n" + about_plugin_txt + "\n<!-- /about_plugin -->" + readme_md_2
    )

    # Write to output file
    success, err_msg = write_file(f"README/{filename}", readme_plugin_md, return_err=True)
    if success:
        output_text(FLAG_SUCCESS)
        output_text(f"<soft>Updated</soft> <reset>/README/{filename}</reset>")
    else:
        output_text(FLAG_ERROR)
        output_error(err_msg, pad=0)


# endregion

############################################################
# region - README/commands.md (GitHub)


# Update the README/commands.md with auto-generated commands
def generate_github_readme_commands_md():
    generate_commands_md("commands.md", for_github=True)


# endregion

############################################################
# region - README/README_pypi.md (PyPi)


# Create the README_pypi.md for PyPI package page
# - - -
# A modified version of the README.md file with links and images
# pointing to the documentation website instead of other readme files.
def generate_pypi_readme_md(filename="README_pypi.md"):
    output_text(f"<h1>Generating <yellow>{filename}</yellow></h1>", pad_top=2)

    # Read README.md file content
    readme_md, err_msg = open_file("README.md", return_err=True)
    if not readme_md:
        output_text(FLAG_ERROR, pad_top=1)
        output_error(err_msg)
        return

    # Remove comments from README content
    readme_md = re.sub(r"<!--.*?-->\n?", "", readme_md, flags=re.DOTALL)

    # Translate image URLs
    readme_md = _translate_image_urls(readme_md)

    # Translate link pointers
    readme_md = _translate_links(readme_md, absolute=True)

    # Insert DO NOT EDIT comment
    readme_md = DO_NOT_EDIT_PYPI + "\n\n" + readme_md

    # Write to output file
    success, err_msg = write_file(f"README/{filename}", readme_md, return_err=True)
    if success:
        output_text(FLAG_SUCCESS)
        output_text(f"<soft>Exported to</soft> <reset>/{filename}</reset>")
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
    base_md_path = "README.md" if filename == "index.md" else f"README/{filename}"
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

    # Remove superfluous linebreaks
    base_md = re.sub(r"\n\s*\n", "\n\n", base_md)

    # Trim space from start and end of file
    base_md = base_md.strip()

    # Translate image URLs
    base_md = _translate_image_urls(base_md)

    # Update link pointers
    base_md = _translate_links(base_md)

    # Update GitHub alerts
    base_md = _translate_alerts(base_md)

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


def _translate_image_urls(readme_md):
    base = "https://raw.githubusercontent.com/acceleratedscience/open-ad-toolkit/main/assets/"

    # Markdown URLs
    # ![some alt text](assets/the_image.png)
    readme_md = re.sub(
        r"!\[(.*?)\]\(assets/(.*?)\)",
        lambda m: f"![{m.group(1)}]({base}{m.group(2)})",
        readme_md,
    )

    # HTML URLs
    # <img src="assets/the_image.png">
    readme_md = re.sub(
        r'<img src="assets/(.*?)"',
        lambda m: f'<img src="{base}{m.group(1)}"',
        readme_md,
    )

    return readme_md


def _translate_links(text, absolute=False):
    """
    When translating a GitHub README file to a just-the-docs markdown file,
    we need to update all internal links so they point to the just-the-docs
    page (xxx.html) instead of the README file (README_xxx.md).

    Similarly, we need to update all internal links for the PyPI README file,
    so they point to the absolute URLs of the documentation website.

    Parameters
    ----------
    text : str
        The text in which to replace the links
    absolute : bool
        If True, the links will be replaced with absolute URLs (for PyPI    )
    """

    absolute_prefix = "https://acceleratedscience.github.io/openad-docs/" if absolute else ""

    # Update template links to the main README.md
    #
    # Input:
    #   [some text]: README.md
    #   [some text]: /README.md
    #   [some text]: ../README.md
    #   [some text]: README.md#foo
    #   [some text]: /README.md#foo
    #   [some text]: ../README.md#foo
    #
    # Relative output (docs website):
    #   [some text]: index.html
    #   [some text]: index.html#foo
    #
    # Absolute output (PyPI):
    #   [some text]: https://acceleratedscience.github.io/openad-docs/
    #   [some text]: https://acceleratedscience.github.io/openad-docs/#foo
    home_link = absolute_prefix if absolute else "index.html"
    text = re.sub(
        r"^\[(.*?)\]: ((\.\./)|(/))?README.md(#.*)?",
        lambda m: f"[{m.group(1)}]: {home_link}{m.group(5) or ''}",
        text,
        flags=re.MULTILINE,
    )

    # Update inline links to the main README.md
    #
    # Input:
    #   [some text](README.md)
    #   [some text](/README.md)
    #   [some text](../README.md)
    #   [some text](README.md#foo)
    #   [some text](/README.md#foo)
    #   [some text](../README.md#foo)
    #
    # Relative output (docs website):
    #   [some text](index.html)
    #   [some text](index.html#foo)
    #
    # Absolute output (PyPI):
    #   [some text](https://acceleratedscience.github.io/openad-docs/)
    #   [some text](https://acceleratedscience.github.io/openad-docs/#foo)
    text = re.sub(
        r"\[(.*?)\]\(((\.\./)|(/))?README.md(#.*)?\)", lambda m: f"[{m.group(1)}]({home_link}{m.group(5) or ''})", text
    )

    # Update template links to secondary README/xxx.md pages
    #
    # Input:
    #   [some text]: abc.md
    #   [some text]: /abc.md
    #   [some text]: README/abc.md
    #   [some text]: /README/abc.md
    #   [some text]: abc.md#foo
    #   [some text]: /abc.md#foo
    #   [some text]: README/abc.md#foo
    #   [some text]: /README/abc.md#foo
    #
    # Relative output (docs website):
    #   [some text]: abc.html
    #   [some text]: abc.html#foo
    #
    # Absolute output (PyPI):
    #   [some text]: https://acceleratedscience.github.io/openad-docs/abc.html
    #   [some text]: https://acceleratedscience.github.io/openad-docs/abc.html#foo
    text = re.sub(
        r"^\[(.*?)\]: (/?README)?/?(.*?).md(#.*)?",
        lambda m: f"[{m.group(1)}]: {absolute_prefix}{m.group(3)}.html{m.group(4) or ''}",
        text,
        flags=re.MULTILINE,
    )

    # Update inline links to secondary README/xxx.md pages
    #
    # Input:
    #   [some text](_abc.md)
    #   [some text](/abc.md)
    #   [some text](README/abc.md)
    #   [some text](/README/abc.md)
    #   [some text](abc.md#foo)
    #   [some text](/abc.md#foo)
    #   [some text](README/abc.md#foo)
    #   [some text](/README/abc.md#foo)
    #
    # Relative output (docs website):
    #   [some text](abc.html)
    #   [some text](abc.html#foo)
    #
    # Absolute output (PyPI):
    #   [some text](https://acceleratedscience.github.io/openad-docs/abc.html)
    #   [some text](https://acceleratedscience.github.io/openad-docs/abc.html#foo)
    text = re.sub(
        r"\[(.*?)\]\((/?README)?/?(.*?).md(#.*)?\)",
        lambda m: f"[{m.group(1)}]({absolute_prefix}{m.group(3)}.html{m.group(4) or ''})",
        text,
    )

    return text


def _translate_alerts(text):
    """
    We use GitHub-specific alert styles in the README files.
    These need to be tranbslated into HTML to render nicely on Just the Docs and PyPI.

    Documentation:
    https://docs.github.com/en/get-started/writing-on-github/getting-started-with-writing-and-formatting-on-github/basic-writing-and-formatting-syntax#alerts
    """

    # Their colors below are the same colors used by GitHub
    color_note = "#0969da"
    color_tip = "#1a7f37"
    color_important = "#8250df"
    color_warning = "#9a6700"
    color_caution = "#d1242f"

    # The icons below are the same SVG icons used in the openad-vue-template repository: https://github.com/acceleratedscience/openad-vue-template
    icn_info = f'<svg class="alert-icon" width="16" height="16" viewBox="0 0 16 16" fill="{color_note}" xmlns="http://www.w3.org/2000/svg"><path d="M8.5 11V7H6.5V8H7.5V11H6V12H10V11H8.5Z"/><path d="M8 4C7.85167 4 7.70666 4.04399 7.58333 4.1264C7.45999 4.20881 7.36386 4.32595 7.30709 4.46299C7.25033 4.60004 7.23547 4.75084 7.26441 4.89632C7.29335 5.04181 7.36478 5.17544 7.46967 5.28033C7.57456 5.38522 7.7082 5.45665 7.85369 5.48559C7.99917 5.51453 8.14997 5.49968 8.28701 5.44291C8.42406 5.38615 8.54119 5.29002 8.6236 5.16668C8.70602 5.04334 8.75 4.89834 8.75 4.75C8.75 4.55109 8.67098 4.36033 8.53033 4.21967C8.38968 4.07902 8.19892 4 8 4Z"/><path d="M8 15C6.61553 15 5.26216 14.5895 4.11101 13.8203C2.95987 13.0511 2.06266 11.9579 1.53285 10.6788C1.00303 9.3997 0.86441 7.99224 1.13451 6.63437C1.4046 5.2765 2.07129 4.02922 3.05026 3.05026C4.02922 2.07129 5.2765 1.4046 6.63437 1.13451C7.99224 0.86441 9.3997 1.00303 10.6788 1.53285C11.9579 2.06266 13.0511 2.95987 13.8203 4.11101C14.5895 5.26216 15 6.61553 15 8C15 9.85652 14.2625 11.637 12.9497 12.9497C11.637 14.2625 9.85652 15 8 15ZM8 2C6.81332 2 5.65328 2.3519 4.66658 3.01119C3.67989 3.67047 2.91085 4.60755 2.45673 5.7039C2.0026 6.80026 1.88378 8.00666 2.11529 9.17054C2.3468 10.3344 2.91825 11.4035 3.75736 12.2426C4.59648 13.0818 5.66558 13.6532 6.82946 13.8847C7.99335 14.1162 9.19975 13.9974 10.2961 13.5433C11.3925 13.0892 12.3295 12.3201 12.9888 11.3334C13.6481 10.3467 14 9.18669 14 8C14 6.4087 13.3679 4.88258 12.2426 3.75736C11.1174 2.63214 9.5913 2 8 2Z"/></svg>'
    icn_bulb = f'<svg class="alert-icon" width="16" height="16" viewBox="0 0 16 16" fill="{color_tip}" xmlns="http://www.w3.org/2000/svg"><path d="M10.5 12H5.49999V13H10.5V12Z"/><path d="M9.49999 14H6.49999V15H9.49999V14Z"/><path d="M7.99999 1C6.67391 1 5.40214 1.52678 4.46446 2.46447C3.52678 3.40215 2.99999 4.67392 2.99999 6C2.96618 6.72667 3.10538 7.45098 3.40614 8.11335C3.7069 8.77572 4.16063 9.35722 4.72999 9.81C5.22999 10.275 5.49999 10.54 5.49999 11H6.49999C6.49999 10.08 5.94499 9.565 5.40499 9.07C4.93767 8.71213 4.56523 8.24514 4.32028 7.70992C4.07534 7.1747 3.96536 6.58759 3.99999 6C3.99999 4.93913 4.42142 3.92172 5.17157 3.17157C5.92171 2.42143 6.93913 2 7.99999 2C9.06086 2 10.0783 2.42143 10.8284 3.17157C11.5786 3.92172 12 4.93913 12 6C12.034 6.58802 11.9233 7.17541 11.6775 7.71067C11.4316 8.24592 11.0582 8.71267 10.59 9.07C10.055 9.57 9.49999 10.07 9.49999 11H10.5C10.5 10.54 10.765 10.275 11.27 9.805C11.839 9.35299 12.2925 8.77235 12.5932 8.11085C12.894 7.44935 13.0334 6.72589 13 6C13 5.34339 12.8707 4.69321 12.6194 4.08658C12.3681 3.47995 11.9998 2.92876 11.5355 2.46447C11.0712 2.00017 10.52 1.63188 9.91341 1.3806C9.30678 1.12933 8.6566 1 7.99999 1Z"/></svg>'
    icn_important = f'<svg class="alert-icon" width="16" height="16" viewBox="0 0 16 16" fill="{color_important}" xmlns="http://www.w3.org/2000/svg"><path d="M8 1C6.61553 1 5.26216 1.41054 4.11101 2.17971C2.95987 2.94888 2.06266 4.04213 1.53285 5.32122C1.00303 6.6003 0.86441 8.00776 1.13451 9.36563C1.4046 10.7235 2.07129 11.9708 3.05026 12.9497C4.02922 13.9287 5.2765 14.5954 6.63437 14.8655C7.99224 15.1356 9.3997 14.997 10.6788 14.4672C11.9579 13.9373 13.0511 13.0401 13.8203 11.889C14.5895 10.7378 15 9.38447 15 8C15 6.14348 14.2625 4.36301 12.9497 3.05025C11.637 1.7375 9.85652 1 8 1ZM8 14C6.81332 14 5.65328 13.6481 4.66658 12.9888C3.67989 12.3295 2.91085 11.3925 2.45673 10.2961C2.0026 9.19974 1.88378 7.99334 2.11529 6.82946C2.3468 5.66557 2.91825 4.59647 3.75736 3.75736C4.59648 2.91824 5.66558 2.3468 6.82946 2.11529C7.99335 1.88378 9.19975 2.0026 10.2961 2.45672C11.3925 2.91085 12.3295 3.67988 12.9888 4.66658C13.6481 5.65327 14 6.81331 14 8C14 9.5913 13.3679 11.1174 12.2426 12.2426C11.1174 13.3679 9.5913 14 8 14Z"/><path d="M8.5 4H7.5V9.5H8.5V4Z"/><path d="M8 11C7.85167 11 7.70666 11.044 7.58333 11.1264C7.45999 11.2088 7.36386 11.3259 7.30709 11.463C7.25033 11.6 7.23547 11.7508 7.26441 11.8963C7.29335 12.0418 7.36478 12.1754 7.46967 12.2803C7.57456 12.3852 7.7082 12.4567 7.85369 12.4856C7.99917 12.5145 8.14997 12.4997 8.28701 12.4429C8.42406 12.3861 8.54119 12.29 8.6236 12.1667C8.70602 12.0433 8.75 11.8983 8.75 11.75C8.75 11.5511 8.67098 11.3603 8.53033 11.2197C8.38968 11.079 8.19892 11 8 11Z"/></svg>'
    icn_warning = f'<svg class="alert-icon" width="16" height="16" viewBox="0 0 16 16" fill="{color_warning}" xmlns="http://www.w3.org/2000/svg"><path d="M8 11.5C7.85167 11.5 7.70666 11.544 7.58333 11.6264C7.45999 11.7088 7.36386 11.826 7.30709 11.963C7.25033 12.1 7.23548 12.2508 7.26441 12.3963C7.29335 12.5418 7.36478 12.6754 7.46967 12.7803C7.57456 12.8852 7.7082 12.9567 7.85369 12.9856C7.99917 13.0145 8.14997 12.9997 8.28702 12.9429C8.42406 12.8862 8.54119 12.79 8.62361 12.6667C8.70602 12.5433 8.75 12.3983 8.75 12.25C8.75 12.0511 8.67099 11.8603 8.53033 11.7197C8.38968 11.579 8.19892 11.5 8 11.5Z"/><path d="M8.5 6.00001H7.5V10.5H8.5V6.00001Z"/><path d="M14.5 15H1.5C1.4141 15 1.32965 14.9779 1.25478 14.9357C1.17992 14.8936 1.11718 14.8329 1.0726 14.7595C1.02802 14.686 1.00311 14.6024 1.00027 14.5165C0.997436 14.4307 1.01677 14.3455 1.0564 14.2693L7.5564 1.76931C7.59862 1.68812 7.66231 1.62008 7.74053 1.5726C7.81875 1.52511 7.9085 1.5 8 1.5C8.09151 1.5 8.18126 1.52511 8.25948 1.5726C8.3377 1.62008 8.40138 1.68812 8.4436 1.76931L14.9436 14.2693C14.9832 14.3455 15.0026 14.4307 14.9997 14.5165C14.9969 14.6024 14.972 14.686 14.9274 14.7595C14.8828 14.8329 14.8201 14.8936 14.7452 14.9357C14.6704 14.9779 14.5859 15 14.5 15ZM2.32535 14H13.6747L13.6757 13.9984L8.001 3.08571H7.999L2.32435 13.9984L2.32535 14Z"/></svg>'
    icn_caution = f'<svg class="alert-icon" width="16" height="16" viewBox="0 0 16 16" fill="{color_caution}" xmlns="http://www.w3.org/2000/svg"><path d="M8.00002 10.5C7.85168 10.5 7.70668 10.544 7.58334 10.6264C7.46 10.7088 7.36388 10.8259 7.30711 10.963C7.25034 11.1 7.23549 11.2508 7.26443 11.3963C7.29337 11.5418 7.3648 11.6754 7.46969 11.7803C7.57458 11.8852 7.70822 11.9567 7.8537 11.9856C7.99919 12.0145 8.14999 11.9997 8.28703 11.9429C8.42408 11.8861 8.54121 11.79 8.62362 11.6667C8.70603 11.5433 8.75002 11.3983 8.75002 11.25C8.75002 11.0511 8.671 10.8603 8.53035 10.7197C8.3897 10.579 8.19893 10.5 8.00002 10.5Z"/><path d="M8.50002 4H7.50002V9H8.50002V4Z"/><path d="M11.5 14.5H4.50002C4.41263 14.5 4.32677 14.4771 4.25099 14.4336C4.1752 14.3901 4.11215 14.3274 4.06812 14.252L0.568119 8.25195C0.523507 8.17548 0.5 8.08853 0.5 8C0.5 7.91147 0.523507 7.82452 0.568119 7.74805L4.06812 1.74805C4.11215 1.67257 4.1752 1.60994 4.25099 1.56642C4.32677 1.5229 4.41263 1.5 4.50002 1.5H11.5C11.5874 1.5 11.6733 1.5229 11.7491 1.56642C11.8248 1.60994 11.8879 1.67257 11.9319 1.74805L15.4319 7.74805C15.4765 7.82452 15.5 7.91147 15.5 8C15.5 8.08853 15.4765 8.17548 15.4319 8.25195L11.9319 14.252C11.8879 14.3274 11.8248 14.3901 11.7491 14.4336C11.6733 14.4771 11.5874 14.5 11.5 14.5ZM4.78712 13.5H11.2129L14.4212 8L11.2129 2.5H4.78712L1.57887 8L4.78712 13.5Z"/></svg>'

    # Replace NOTE blocks
    note_pattern = r"> \[!NOTE\]((?:\n> .*)*)"
    note_matches = re.findall(note_pattern, text, flags=re.MULTILINE)
    for match in note_matches:
        match = re.sub(r"^\n> ", "", match).replace("\n> ", "<br>")
        text = re.sub(
            note_pattern,
            f"> <div class='alert-icn-wrap' style='color:{color_note}'>{icn_info} NOTE</div><span style='color: {color_note}'>{match}</span>",
            text,
            flags=re.MULTILINE,
            count=1,
        )

    # Replace TIP blocks
    tip_pattern = r"> \[!TIP\]((?:\n> .*)*)"
    tip_matches = re.findall(tip_pattern, text, flags=re.MULTILINE)
    for match in tip_matches:
        match = re.sub(r"^\n> ", "", match).replace("\n> ", "<br>")
        text = re.sub(
            tip_pattern,
            f"> <div class='alert-icn-wrap' style='color:{color_tip}'>{icn_bulb} TIP</div><span style='color: {color_tip}'>{match}</span>",
            text,
            flags=re.MULTILINE,
            count=1,
        )

    # Replace IMPORTANT blocks
    important_pattern = r"> \[!IMPORTANT\]((?:\n> .*)*)"
    important_matches = re.findall(important_pattern, text, flags=re.MULTILINE)
    for match in important_matches:
        match = re.sub(r"^\n> ", "", match).replace("\n> ", "<br>")
        text = re.sub(
            important_pattern,
            f"> <div class='alert-icn-wrap' style='color:{color_important}'>{icn_important} IMPORTANT</div><span style='color: {color_important}'>{match}</span>",
            text,
            flags=re.MULTILINE,
            count=1,
        )

    # Replace WARNING blocks
    warning_pattern = r"> \[!WARNING\]((?:\n> .*)*)"
    warning_matches = re.findall(warning_pattern, text, flags=re.MULTILINE)
    for match in warning_matches:
        match = re.sub(r"^\n> ", "", match).replace("\n> ", "<br>")
        text = re.sub(
            warning_pattern,
            f"> <div class='alert-icn-wrap' style='color:{color_warning}'>{icn_warning} WARNING</div><span style='color: {color_warning}'>{match}</span>",
            text,
            flags=re.MULTILINE,
            count=1,
        )

    # Replace CAUTION blocks
    caution_pattern = r"> \[!CAUTION\]((?:\n> .*)*)"
    caution_matches = re.findall(caution_pattern, text, flags=re.MULTILINE)
    for match in caution_matches:
        match = re.sub(r"^\n> ", "", match).replace("\n> ", "<br>")
        text = re.sub(
            caution_pattern,
            f"> <div class='alert-icn-wrap' style='color:{color_caution}'>{icn_caution} CAUTION</div><span style='color: {color_caution}'>{match}</span>",
            text,
            flags=re.MULTILINE,
            count=1,
        )

    return text


# endregion

############################################################
# region - base-concepts.md


# Generate the base-concepts.md file for the documentation website.
def generate_base_concepts_md(filename="base-concepts.md"):
    output_text(f"<h1>Generating <yellow>{filename}</yellow></h1>", pad_top=2)

    # Read base-concepts.md input file content
    base_concepts_md, err_msg = open_file("docs/input/base-concepts.md", return_err=True)
    if not base_concepts_md:
        output_text(FLAG_ERROR, pad_top=1)
        output_error(err_msg)
        return

    # Read about_workspace.txt source file content
    about_workspace, err_msg = open_file("openad/docs_src/about_workspace.txt", return_err=True)
    if not about_workspace:
        output_text(FLAG_ERROR, pad_top=1)
        output_error(err_msg)
        return

    # Read about_mws.txt source file content (molecule working set)
    about_mws, err_msg = open_file("openad/docs_src/about_mws.txt", return_err=True)
    if not about_mws:
        output_text(FLAG_ERROR, pad_top=1)
        output_error(err_msg)
        return

    # Read about_plugin.txt source file content
    about_plugin, err_msg = open_file("openad/docs_src/about_plugin.txt", return_err=True)
    if not about_plugin:
        output_text(FLAG_ERROR, pad_top=1)
        output_error(err_msg)
        return

    # Read about_context.txt source file content
    about_context, err_msg = open_file("openad/docs_src/about_context.txt", return_err=True)
    if not about_context:
        output_text(FLAG_ERROR, pad_top=1)
        output_error(err_msg)
        return

    # Read about_run.txt source file content
    about_run, err_msg = open_file("openad/docs_src/about_run.txt", return_err=True)
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
    success, err_msg = write_file(f"docs/output/markdown/{filename}", base_concepts_md, return_err=True)
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
def generate_commands_md(filename="commands.md", for_github=False):
    output_text(f"<h1>Generating <yellow>{filename}</yellow> from help</h1>", pad_top=2)

    toc = []  # Table of content
    md_output = []  # Markdown

    # Parse main commands
    space = "<br>\n\n" if for_github else "<br><br>\n\n"
    md_output.append(f"{space}## Main Commands\n")
    toc.append(_toc_link("Main Commands"))
    cmds = cmd_pointer.current_help.help_current
    cmds_organized = organize_commands(cmds)
    if for_github:
        _compile_section_github(md_output, toc, cmds_organized)
    else:
        _compile_section(md_output, toc, cmds_organized)

    # Parse tookit commands
    for toolkit_name in _all_toolkits:
        space = "<br>\n\n" if for_github else "<br><br>\n\n"
        md_output.append(f"{space}## {toolkit_name}\n\n")
        toc.append(_toc_link(toolkit_name))
        success, toolkit = load_toolkit(toolkit_name, from_repo=True)
        if success:
            toolkit_cmds = toolkit.methods_help
            toolkit_cmds_organized = organize_commands(toolkit_cmds)
            if for_github:
                _compile_section_github(md_output, toc, toolkit_cmds_organized)
            else:
                _compile_section(md_output, toc, toolkit_cmds_organized)

    # Compile table of contents
    toc = "<br>\n\n## Table of Contents\n" + "\n".join(toc) + "\n"

    # Compile commands
    md_output = "\n".join(md_output)

    # Read commands.md input content
    commands_md, err_msg = open_file("docs/input/commands.md", return_err=True)
    if not commands_md:
        output_text(FLAG_ERROR, pad_top=1)
        output_error(err_msg)
        return

    # Update GitHub alerts
    if not for_github:
        commands_md = _translate_alerts(commands_md)

    # Replace the just-the-docs header with a back link when generating for GitHub
    if for_github:
        commands_md = re.sub(r"^---.+---", "<sub>[&larr; BACK](../#openad)</sub>", commands_md, flags=re.DOTALL)

    # Insert DO NOT EDIT comment
    commands_md = re.sub(r"{{DO_NOT_EDIT}}", DO_NOT_EDIT, commands_md, flags=re.DOTALL)

    # Insert table of contents
    commands_md = re.sub(r"{{TOC}}", toc, commands_md, flags=re.DOTALL)

    # Insert commands
    commands_md = re.sub(r"{{COMMANDS}}", md_output, commands_md, flags=re.DOTALL)

    # Write to file
    if for_github:
        success, err_msg = write_file(f"README/{filename}", commands_md, return_err=True)
    else:
        success, err_msg = write_file(f"docs/output/markdown/{filename}", commands_md, return_err=True)
    if success:
        output_text(FLAG_SUCCESS)
        if for_github:
            output_text(f"<soft>Exported to</soft> <reset>/README/{filename}</reset>")
        else:
            output_text(f"<soft>Exported to</soft> <reset>/docs/output/markdown/{filename}</reset>")
    else:
        output_text(FLAG_ERROR)
        output_error(err_msg, pad=0)


# Compile all commands of a single section.
def _compile_section(output, toc, cmds_organized):
    for category in cmds_organized:
        output.append(f"### {category}\n")
        toc.append(_toc_link(category, 1))
        for cmd_str, cmd_description in cmds_organized[category]:
            cmd_output = "\n".join(
                [
                    '<details markdown="block" class="cmd-wrap">',
                    '<summary markdown="block">',
                    f"`{cmd_str.strip()}`{{: .cmd }}",
                    "</summary>",
                    _parse_description(cmd_description),
                    "</details>\n",
                ]
            )
            output.append(cmd_output)


def _compile_section_github(output, toc, cmds_organized):
    for i, category in enumerate(cmds_organized):
        space = "" if i == 0 else "<br>\n\n"
        output.append(f"{space}### {category}\n")
        toc.append(_toc_link(category, 1))
        for cmd_str, cmd_description in cmds_organized[category]:
            # Add `> ` in front of every line so the description shows up as a note block
            cmd_description = (
                "<br>\n\n" + re.sub(r"^", "> ", _parse_description(cmd_description), flags=re.MULTILINE) + "\n"
            )
            cmd_output = "\n".join(
                [
                    '<details markdown="block" class="cmd-wrap">',
                    f'<summary markdown="block"><code>{cmd_str.strip()}</code></summary>',
                    cmd_description,
                    "</details>\n",
                ]
            )
            output.append(cmd_output)


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
def generate_commands_csv(filename="commands.csv", delimiter=";"):
    output_text("<h1>Generating <yellow>commands.csv</yellow> from help</h1>", pad_top=2)
    output = [["Command", "Category"]]

    # Parse main commands
    cmds_main = cmd_pointer.current_help.help_current
    cmds_organized = organize_commands(cmds_main)

    # Parse tookit commands
    for toolkit_name in _all_toolkits:
        success, toolkit = load_toolkit(toolkit_name, from_repo=True)
        if success:
            toolkit_cmds = toolkit.methods_help
            toolkit_cmds_organized = organize_commands(toolkit_cmds)
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
    success, err_msg = write_file(f"docs/output/csv/{filename}", output_str, return_err=True)
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
def generate_llm_description_txt(filename="llm_description.txt"):
    output_text("<h1>Updating commands in <yellow>llm_description.txt</yellow> for all toolkits</h1>", pad_top=2)

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
        toolkit_cmds_organized = organize_commands(toolkit_cmds)
        output = _compile_commands(toolkit_cmds_organized)

        # Load llm_description.txt
        file_path = f"openad/user_toolkits/{toolkit_name}/{filename}"
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
                f"<soft>Updated in</soft> <reset>/docs/openad/user_toolkits/{toolkit_name}/{filename}</reset>",
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

    # Update existing README files
    output_text("<magenta>Updating existing README files</magenta>", pad_top=4)
    update_github_readme_md()
    update_github_readme_plugins_md()

    # Generate README files
    output_text("<magenta>Generating README files</magenta>", pad_top=4)
    generate_github_readme_commands_md()
    generate_pypi_readme_md()  # For PyPI (links pointing to docs)

    # Turn README files into pages for the documentation website
    output_text("<magenta>Translating README files to doc website pages</magenta>", pad_top=4)
    render_docs_pages()

    # Generate additional bespoke pages for documentation website
    output_text("<magenta>Generate additional doc website pages</magenta>", pad_top=4)
    generate_base_concepts_md()
    generate_commands_md()

    # Render additional files
    output_text("<magenta>Generate additional files</magenta>", pad_top=4)
    generate_commands_csv()
    generate_llm_description_txt()

    # For testing
    # _render_docs_page("index.md")

    # Move all generated markdown files to the documentation repo
    docs = []
    for filename in os.listdir(f"docs/output/markdown"):
        docs.append(filename)
    copy_docs(docs)
