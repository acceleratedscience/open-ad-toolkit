import openad.app.main  # we import entire main library to help retain state
from pandas import DataFrame
from pandas.io.formats.style import Styler
import atexit

context_cache = {"workspace": None, "toolkit": None}


class openad_api:
    """Magic Command Class"""

    def request(self, command):
        """Invokes the Magic command interface for OpenAD and ensure dataFrame Data is of type data"""
        api_variable = {}
        openad.app.main.GLOBAL_SETTINGS["display"] = "api"
        command_list = command.split()
        x = len(command_list)
        i = 1
        if x > 1:
            while i < x:
                if command_list[i - 1].upper() == "DATAFRAME":
                    try:
                        df = eval(command_list[i])  # pylint: disable=eval-used #only way to execute
                        if isinstance(df, DataFrame):
                            api_variable[command_list[i]] = df
                    except:  # pylint: disable=bare-except # We do not care what fails
                        pass
                i += 1

        result = openad.app.main.api_remote(command, context_cache, api_variable)

        if isinstance(result, Styler):
            result = result.data
        return result


def strip_leading_blanks(input):
    temp = input.split("\n")
    output = ""
    for x in temp:
        while str(x).startswith("   "):
            X = str(x).replace("   ", "  ")
        output = output + x + "\n"
    return output


def cleanup():
    print("killing magic")
    openad.app.main.MAGIC_PROMPT.do_exit("exit magic")


atexit.register(cleanup)
