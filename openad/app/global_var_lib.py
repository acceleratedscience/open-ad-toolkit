"""
Central place for commonly used global variables.
"""

import os
import shutil
from openad.helpers.general import get_toolkits
from openad.app.memory import Memory

# client metadata directory
_meta_dir = os.path.expanduser("~/.openad")
_meta_dir_toolkits = os.path.expanduser("~/.openad/toolkits")
_meta_registry = os.path.expanduser("~/.openad/registry.pkl")
_meta_login_registry = os.path.expanduser("~/.openad/login_registry.pkl")
_meta_registry_session = os.path.expanduser("~/.openad/sessions/registry.pkl")
_meta_workspaces = os.path.expanduser("~/.openad/workspaces")
_meta_registry_settings = {
    "workspace": "DEFAULT",
    "context": None,
    "workspaces": ["DEFAULT"],
    "experiments": [],
    "paths": dict(),
    "descriptions": {"DEFAULT": "This is your default workspace."},
    "toolkits": [],
    "env_vars": {},
}
_meta_login_registry_settings = {
    "toolkits": [],
    "toolkits_details": [],
    "toolkits_api": [],
    "client": [],
    "expiry": [],
    "session_vars": [],
}

# Other
_all_toolkits = get_toolkits()
_date_format = "%a %b %d, %G - %R"
_repo_dir = os.path.dirname(os.path.abspath(__file__))


GLOBAL_SETTINGS = {
    # Dictates where our output will be displayed:
    # - terminal: set in main.py -> cmd_line()
    # - notebook: set in main.py -> api_remote()
    # - api: not yet used
    # - web: not yet used
    "display": None,
    #
    # Setting verbose to false will hide spinners
    # and some other non-essential output
    # Todo: this is diplocate of "display" is "notebook" (to be confirmed)
    "VERBOSE": True,
    #
    # The width in characters we limit the printed output to.
    "max_print_width": 150,
    #
    # Flag that triggers the grammar to be refreshed
    # after we add new molecule properties.
    "grammar_refresh": False,
    #
    # "MODEL_SERVICES": None, # unused, trash?
}
MEMORY = Memory()
