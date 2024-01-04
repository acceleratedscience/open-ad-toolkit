#!/usr/local/opt/python@3.9/bin/python3.9
import os
import re
import imp
import glob
import json
import logging

import openad.core.grammar as grammar  # Not using "from" to avoid circular import.
from openad.helpers.output import output_error
from openad.helpers.output_msgs import msg
from openad.helpers.general import load_module_from_path

# Globals
from openad.app.global_var_lib import _meta_dir_toolkits
from openad.app.global_var_lib import _meta_workspaces


# Creates a toolkit object that stores all data about the toolkit.
# - - -
# This function is called to build statements for toolkits or in future other adaptors.
# It transalates a JSON file's contents and creates a statement from it.
# it will be enhanced in future to enable more complex statements as wel work out verbage.
class Toolkit:
    def __init__(self, name) -> None:
        self.toolkit_name = name
        self.toolkit_description = None
        self.methods = []
        self.methods_grammar = []
        self.methods_execute = []
        self.methods_help = []
        self.methods_dict = []
        self.methods_library = []


# Load all toolkit statments.
def load_toolkit(toolkit_name, from_repo=False, for_training=False):
    """
    Parameters
    ----------
    from_repo
        Load the toolkit from the repo instead of the user folder.
        This is used to generate commands.md, where we need to fetch all
        commands from all toolkits, including those that are not installed.

    for_training
        Used by grammar.py to generate the training data.

    """
    # from_repo = True  # TEMP FOR TESTING! DELETE THIS LINE!!

    the_toolkit = Toolkit(toolkit_name)
    if from_repo:
        openad_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
        source = openad_dir + "/user_toolkits"
    else:
        source = _meta_dir_toolkits

    # Load toolkit description snippets.
    snippetsModule = load_module_from_path("snippets", f"{source}/{toolkit_name}/_snippets.py")
    snippets = snippetsModule.snippets if snippetsModule else None

    for i in glob.glob(f"{source}/{toolkit_name}/**/fn_*.json", recursive=True):
        fn_file = open(i, "r", encoding="utf-8")
        x = json.load(fn_file)

        # Load description from separate file if it is external.
        if x["help"]["description"] == "":
            try:
                txt_file = open(i.replace(".json", ".txt"), "r")
                x["help"]["description"] = txt_file.read()
            except BaseException:
                x["help"]["description"] = "Failed to load description"

        # Replace snippet tags with snippet content.
        # - - -
        # We centralize repeating text in one place per toolkit (_snippets.py)
        # which is then referenced in a function's description by tags.
        # For example "lorem ipsum {{FOO_BAR}}" -> "lorem ipsum foo bar"
        if snippets:
            x["help"]["description"] = re.sub(
                r"\{\{(\w+)\}\}",
                lambda match: str(snippets.get(match.group(1), f"[[missing snippet: {match.group(1)}]]").strip()),
                x["help"]["description"],
            )

        grammar.statement_builder(the_toolkit, x)

    # Load toolkit LLM training text.
    if for_training:
        try:
            with open(
                _meta_dir_toolkits + "/" + toolkit_name + "/llm_description.txt", "r", encoding="utf-8"
            ) as toolkit_file:
                the_toolkit.toolkit_description = toolkit_file.read()
                toolkit_file.close()
        except Exception:
            # If unable to load, move on.
            the_toolkit.toolkit_description = None

    return True, the_toolkit


# Load toolkit description file.
def load_toolkit_description(cmd_pointer, toolkit_name):
    try:
        # We load description from our repo folder instead of the user folder, because
        # we need to be able to access the description before the toolkit is installed.
        # (list all toolkits)
        # desc_file = open(_meta_dir_toolkits + "/" + toolkit_name + "/oneline_desc.txt", "r")
        desc_file = open(
            os.path.dirname(os.path.abspath(__file__)) + "/../user_toolkits/" + toolkit_name + "/oneline_desc.txt", "r"
        )
        return str(desc_file.readline())
    except BaseException as err:
        return output_error(msg("err_workspace_description"), return_val=True)


# Enact a method call to a toolkit library.
def execute_tookit(cmd_pointer, parser):
    name = parser.getName().replace("toolkit_exec_", "")
    index = cmd_pointer.toolkit_current.methods.index(name)

    logging.basicConfig(
        level=logging.INFO, filename=_meta_workspaces + "/workspace.log", format="%(asctime)s - %(message)s"
    )

    for i in glob.glob(
        _meta_dir_toolkits
        + "/"
        + cmd_pointer.toolkit_current.toolkit_name
        + "/**/"
        + cmd_pointer.toolkit_current.methods_library[index]
        + ".py",
        recursive=True,
    ):
        exec_link = load_src(cmd_pointer.toolkit_current.methods_library[index], i)
        func = getattr(exec_link, cmd_pointer.toolkit_current.methods_execute[index])
        test_typing(parser, cmd_pointer.toolkit_current.methods_dict[index])

        logging.info(
            "Executing "
            + cmd_pointer.toolkit_current.methods[index]
            + " From "
            + cmd_pointer.toolkit_current.methods_library[index]
        )
        try:
            x = func(parser.as_dict(), cmd_pointer)
            logging.info(
                "Executed "
                + cmd_pointer.toolkit_current.methods[index]
                + " From "
                + cmd_pointer.toolkit_current.methods_library[index]
            )
            return x
        except BaseException as err:
            print(err)
            print("Function call Failed: Likely Server Unresponsive check above error")
            logging.info(
                "Failure "
                + cmd_pointer.toolkit_current.methods[index]
                + " From "
                + cmd_pointer.toolkit_current.methods_library[index]
            )
            return False


# Fetch the load path for a toolkit library.
def load_src(name, fpath):
    return imp.load_source(name, os.path.join(os.path.dirname(__file__), fpath))


# Extra type checking on top of pyparsing, specifically for float types.
# We can add more complex types here if needed.
def test_integer(possible_int) -> bool:
    try:
        int(possible_int) == possible_int
    except:
        return False
    return True


def test_typing(parser, definition):
    results = parser.as_dict()
    for x in results:
        try:
            result = results[x]["val"]
        except BaseException:
            continue
        if isinstance(result, list):
            continue
        if isinstance(result, dict):
            continue
        if "fixed_parameters" in definition:
            if x in definition["fixed_parameters"]:
                if definition["fixed_parameters"][x] == "dec":
                    if result.isnumeric() == False:
                        print("invalid value for paramater " + x)
                        return False
                    else:
                        continue
                elif definition["fixed_parameters"][x] == "int":
                    try:
                        y = int(result)
                        continue
                    except BaseException:
                        print("invalid integer for paramater " + x)
                        return False
                elif definition["fixed_parameters"][x] == "str":
                    if len(result.strip()) == 0:
                        print("invalid string for paramater " + x)
                    else:
                        continue
        # elif  x in definition['variable_parameters']:

        #    if definition['variable_parameters'][x] =='dec':
        #        if result.isnumeric() == False:
        #            print('invalid value for paramater '+x)
        #            return False
        #        else:
        #            continue
        #    elif definition['variable_parameters'][x] =='dec':
        #        try:
        #            y=int(result)
        #            continue
        #        except:
        #            print('invalid integer for paramater '+x)
        #            return False
        #    elif definition['variable_parameters'][x] =='str':
        #            if len(str(result).strip())==0:
        #                print('invalid string for paramater '+x)
        #            else:
        #                continue
