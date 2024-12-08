import importlib
from pandas import DataFrame
from pandas.io.formats.style import Styler
from copy import deepcopy
import re
from openad.app.global_var_lib import _all_toolkits
from openad.toolkit.toolkit_main import load_toolkit
from openad.plugins.style_parser import tags_to_markdown
from openad.helpers.output import output_error, output_text, output_success
from openad.core.help import organize_commands


class OpenadAPI:
    """API class for OpenAD"""

    main_app = None
    module_name = "openad.app.main"
    name = None
    context_cache = deepcopy({"workspace": None, "toolkit": None})

    def __init__(self, name="No Name"):
        # import openad.app.main as main_app

        self.main_app = self._load_main()
        self.main_app.GLOBAL_SETTINGS["VERBOSE"] = False
        self.name = name

    def _load_main(self):
        spec = importlib.util.find_spec(self.module_name)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        return module

    def request(self, command, Vars=None, **kwargs):
        """Invokes the Magic command interface for OpenAD and ensure dataFrame Data is of type data"""
        api_variable = {}
        self.main_app.GLOBAL_SETTINGS["display"] = "api"
        command_list = command.split()
        x = len(command_list)
        i = 1
        if x > 1:
            while i < x:
                if command_list[i - 1].upper() == "DATAFRAME":
                    try:
                        df = kwargs[command_list[i]]  # pylint: disable=eval-used #only way to execute
                        if isinstance(df, DataFrame):
                            api_variable[command_list[i]] = df
                    except:  # pylint: disable=bare-except # We do not care what fails
                        pass
                i += 1

        result = self.main_app.api_remote(command, self.context_cache, api_variable)

        if isinstance(result, Styler):
            result = result.data

        return result

    def help_as_markdown(self, command):
        x = self.main_app.RUNCMD().do_help(command, display_info=False, jup_return_format=True)
        return x

    def __del__(self):
        return

    def help_dump(self):
        """dumps the help text in markup"""

        # return render_commands_csv()

        output_text("<h1>Generating <yellow>commands.md</yellow> from help</h1>", pad_top=4)

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
            "To update it, see openad/docs/generate_docs.py",
        )
        comment = "\n".join(comment)
        output.append(f"<!--\n\n{comment}\n\n-->" + "\n")

        # Parse main commands
        output.append(f"## OpenAD\n")
        toc.append(_toc_link("OpenAD"))
        cmds = self.main_app.RUNCMD().current_help.help_current
        cmds_organized = organize_commands(cmds)
        _compile_section(output, toc, cmds_organized)

        # Parse tookit commands
        for toolkit_name in _all_toolkits:
            output.append(f"## {toolkit_name}\n\n")
            toc.append(_toc_link(toolkit_name))
            success, toolkit = load_toolkit(toolkit_name, from_repo=True)
            if success:
                toolkit_cmds = toolkit.methods_help
                toolkit_cmds_organized = organize_commands(toolkit_cmds)
                _compile_section(output, toc, toolkit_cmds_organized)

        # Write output to file to this python file's parent folder
        toc = "### Table of Contents\n" + "\n".join(toc) + "\n"
        output = output[:2] + [toc] + output[2:]
        output = "\n".join(output)
        return output

    def strip_leading_blanks(self, input):
        temp = input.split("\n")
        output = ""
        for x in temp:
            while str(x).startswith("   "):
                X = str(x).replace("   ", "  ")
            output = output + x + "\n"
        return output


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


if __name__ == "__main__":
    myclass = OpenadAPI("class1")
