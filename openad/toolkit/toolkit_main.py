#!/usr/local/opt/python@3.9/bin/python3.9
import os
import imp
import glob
import json
import logging
from openad.core.grammar import statement_builder
from openad.helpers.output import msg, output_error

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
def load_toolkit(toolkit_name):
    the_toolkit = Toolkit(toolkit_name)

    for i in glob.glob(_meta_dir_toolkits + "/" + toolkit_name + "/**/func_*.json", recursive=True):
        func_file = open(i, "r")
        x = json.load(func_file)
        statement_builder(the_toolkit, x)

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
        return output_error(msg("err_workspace_description"), cmd_pointer, return_val=True)


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
