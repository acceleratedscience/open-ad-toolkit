# This file has moved to the openad_tools repo.
# All imports should be changed to openad_tools.helpers.
from openad_tools.helpers import *


import os


# To be deleted when we remove toolkits.
# Return list of available toolkit names.
def get_toolkits():
    folder_path = os.path.dirname(os.path.abspath(__file__)) + "/../user_toolkits"
    ignore_dirs = ["__pycache__", "readme"]
    toolkit_names = [
        name.upper()
        for name in os.listdir(folder_path)
        if os.path.isdir(os.path.join(folder_path, name)) and name not in ignore_dirs
    ]
    return toolkit_names
