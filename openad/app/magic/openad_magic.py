import os
import sys
import pandas
import atexit

from openad.app.main import GLOBAL_SETTINGS

# required for Magic Template
from IPython.display import Markdown
from IPython.display import display
from IPython.core.magic import (
    Magics,
    magics_class,
    line_magic,
    cell_magic,
    line_cell_magic,
    needs_local_scope,
)  # pylint: disable=import-error
from IPython.core.interactiveshell import InteractiveShell  # pylint: disable=import-error
import openad.app.main  # we import entire main library to help retain state
from openad.helpers.output import output_table, output_text

InteractiveShell.ast_node_interactivity = "all"
from pandas import DataFrame
from pandas.io.formats.style import Styler

sys.path.insert(0, "../")
os.sys.path.append(os.path.dirname(os.path.abspath("./")))
module_path = os.path.abspath(os.path.join(".."))


handle_cache = {
    "toolkits": [],
    "toolkits_details": [],
    "toolkits_api": [],
    "client": [],
    "expiry": [],
    "session_vars": [],
}
context_cache = {"workspace": None, "toolkit": None}


@magics_class
class AD(Magics):
    """Magic Command Class"""

    @needs_local_scope
    @line_cell_magic
    def openad(self, line, cell=None, local_ns=None):
        """Invokes the Magic command interface for OpenAD"""
        api_variable = {}
        GLOBAL_SETTINGS["display"] = "notebook"
        line_list = line.split()
        x = len(line_list)
        i = 1
        if x > 1:
            while i < x:
                if line_list[i - 1].upper() == "DATAFRAME":
                    try:
                        df = eval(line_list[i])  # pylint: disable=eval-used #only way to execute
                        if isinstance(df, pandas.DataFrame):
                            api_variable[line_list[i]] = df
                    except:  # pylint: disable=bare-except # We do not care what fails
                        pass
                i += 1
        result = openad.app.main.api_remote(line, context_cache, api_variable)

        if isinstance(result, DataFrame):
            result = output_table(result, return_val=True)
        elif isinstance(result, str):
            result = strip_leading_blanks(result)
            result = result.replace("<br>", "\n")

        return result

    @needs_local_scope
    @line_cell_magic
    def openadd(self, line, cell=None, local_ns=None):
        """Invokes the Magic command interface for OpenAD and ensure dataFrame Data is of type data"""
        api_variable = {}
        GLOBAL_SETTINGS["display"] = "api"
        line_list = line.split()
        x = len(line_list)
        i = 1
        if x > 1:
            while i < x:
                if line_list[i - 1].upper() == "DATAFRAME":
                    try:
                        df = eval(line_list[i])  # pylint: disable=eval-used #only way to execute
                        if isinstance(df, pandas.DataFrame):
                            api_variable[line_list[i]] = df
                    except:  # pylint: disable=bare-except # We do not care what fails
                        pass
                i += 1
        result = openad.app.main.api_remote(line, context_cache, api_variable)

        # MAJOR-RELEASE-TODO: data function should never display
        if isinstance(result, Styler):
            result = result.data
        # if isinstance(result, str):
        #    display(Markdown(result))
        return result


def strip_leading_blanks(input):
    temp = input.split("\n")
    output = ""
    for x in temp:
        while str(x).startswith("   "):
            X = str(x).replace("   ", "  ")
        output = output + x + "\n"
    return output


ip = get_ipython()  # pylint: disable=undefined-variable
ip.register_magics(AD)


def cleanup():
    print("killing magic")
    openad.app.main.MAGIC_PROMPT.do_exit("exit magic")


atexit.register(cleanup)
