"""Used for logging into toolkits"""

# python 3.8
import pickle
import os
import imp  # TODO: future improvement to upgrade from imp
from openad.app.global_var_lib import _meta_login_registry
from openad.app.global_var_lib import _meta_login_registry_settings
from openad.helpers.output import output_success, output_error
from openad.helpers.output_msgs import msg
from openad.app.global_var_lib import _meta_dir_toolkits

# Initialise  Login Pickle


def initialise_toolkit_login():
    """Initialises the settings file for login, currently this code
    is a holding place for a persist to disk solution for storing login handles"""
    try:
        with open(_meta_login_registry, "wb") as handle:
            pickle.dump(_meta_login_registry_settings, handle)
            handle.close()
        output_success(msg("success_login_init"), return_val=False)
        return True
    except Exception as err:
        output_error(msg("error_login_init", err), return_val=False)
        return False


# loads the user's registry data
def load_login_registry():
    """Loads connection handles from disk, currently this code is
    a holding place for a persist to disk solution for storing login handles"""
    try:
        with open(_meta_login_registry, "rb") as handle:
            login_registry = pickle.load(handle.read())
            handle.close()
    except Exception:
        login_registry = _meta_login_registry_settings.copy()
    return login_registry


# writes the users registry data, this is dummy function until we work through handl caching on disk...
def write_login_registry(login_registry: dict):
    """writes connection handles to disk, currently this code is
    a holding place for a persist to disk solution for storing login handles"""
    return True


# writes the users registry data
def load_src(name, fpath):
    """loads the source library"""
    return imp.load_source(name, os.path.join(os.path.dirname(__file__), fpath))


def load_login_api(cmd_pointer, toolkit_name, reset=False):
    """loads the login API"""
    suppress = False
    if toolkit_name.upper() in cmd_pointer.settings["toolkits"]:
        suppress = True
    if os.path.isfile(f"{_meta_dir_toolkits}/{toolkit_name}/login.py"):
        exec_link = load_src("login", f"{_meta_dir_toolkits}/{toolkit_name}/login.py")
        if reset is True:
            func = getattr(exec_link, "reset")
            func(cmd_pointer)
        func = getattr(exec_link, "login")
        try:
            login_success, expiry_datetime = func(cmd_pointer)
            if expiry_datetime is None:
                expiry_datetime = "No Expiry"
            if login_success is True:
                if not suppress:
                    output_success(msg("success_login", toolkit_name, expiry_datetime), return_val=False)
                return login_success, expiry_datetime
            else:
                if not suppress:
                    output_error(msg("err_login", toolkit_name), return_val=False)
                return False, None

        except Exception as err:
            output_error(msg("err_login", toolkit_name, err), return_val=False)
            return False, None
    else:
        # If no login file is present, assume no login is required.
        return True, None
