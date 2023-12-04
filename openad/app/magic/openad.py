import os
import sys
import pandas

# required for Magic Template
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

InteractiveShell.ast_node_interactivity = "all"


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

        return openad.app.main.api_remote(line, context_cache, api_variable)


ip = get_ipython()  # pylint: disable=undefined-variable
ip.register_magics(AD)
