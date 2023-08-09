# Where global Variables are defined

import os
from helpers.general import get_toolkits

# client metadata directory
_meta_dir = os.path.expanduser("~/.adccl")
_meta_dir_toolkits = os.path.expanduser("~/.adccl/toolkits")
_meta_registry = os.path.expanduser("~/.adccl/registry.pkl")
_meta_login_registry = os.path.expanduser("~/.adccl/login_registry.pkl")
_meta_registry_session = os.path.expanduser("~/.adccl/sessions/registry.pkl")
_meta_workspaces = os.path.expanduser("~/.adccl/workspaces")

_meta_registry_settings = {
    'workspace': 'DEFAULT',
    'context': None,
    'workspaces': ['DEFAULT'],
    'experiments': [],
    'paths': dict(),
    'descriptions': {
        'DEFAULT': 'This is your default workspace.'},
    'toolkits': [],
    'env_vars':{}}
_meta_login_registry_settings = {
    'toolkits': [],
    'toolkits_details': [],
    'toolkits_api': [],
    'client': [],
    'expiry': [],
    'session_vars':[]
}

# Other
_all_toolkits = get_toolkits()
_date_format = '%a %b %d, %G - %R'
#_repo_dir = os.getcwd()  # @Phil the only way I was able to get the repo directory. Feels a bit hacky but it works.
_repo_dir = os.path.dirname(os.path.abspath(__file__))