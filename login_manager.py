# python 3.8
from global_var_lib import _meta_login_registry as _meta_login_registry
from global_var_lib import _meta_login_registry_settings as _meta_login_registry_settings

# Initialise  Login Pickle
import pickle

# The idea here is to give some feedback while the user
# waits a second or so for the auth to complete. But we
# can't properly do this until we centralize login messages.
#
# This line is supposed to be erased as soon as the error/success
# messages returns, using remove_lines(1)
# print(parse_tags('\n<soft>Logging you in...</soft>'))


def initialise_toolkit_login():
    try:
        with open(_meta_login_registry, 'wb') as handle:
            settings = pickle.dump(_meta_login_registry_settings, handle)
            handle.close()
        print("Login registry initiaized.")
        return True
    except BaseException:
        return False


# loads rhe users registry data
def load_login_registry():
    try:
        with open(_meta_login_registry, 'rb') as handle:
            login_registry = pickle.load(handle.read())
            handle.close()
    except BaseException:
        login_registry = _meta_login_registry_settings.copy()
    return login_registry


# writes the users registry data
def write_login_registry(login_registry: dict):
    return True
    print('writing Registry')
    print(login_registry)
    import collections
    try:
        print(login_registry)
        # flush_dict= login_registry.copy()
        with open(_meta_login_registry, 'wb') as handle:
            settings = pickle.dump(login_registry, handle, protocol=None, fix_imports=True, buffer_callback=None)
            handle.close()
        print('wrote Registry')
    except BaseException as err:
        print(err)
        return False
    return True


# writes the users registry data
def load_src(name, fpath):
    import os
    import imp
    return imp.load_source(name, os.path.join(os.path.dirname(__file__), fpath))


def load_login_api(cmd_pointer, toolkit_name,reset=False):
    import os
    from global_var_lib import _meta_dir_toolkits as _meta_dir_toolkits
    suppress = False
    if toolkit_name.upper() in cmd_pointer.settings['toolkits']:
        suppress = True
    if os.path.isfile(_meta_dir_toolkits + "/" + toolkit_name + "/login.py"):
        exec_link = load_src("login", _meta_dir_toolkits + "/" + toolkit_name + "/login.py")
        if reset==True:
            func = getattr(exec_link, "reset")
            func(cmd_pointer)
        func = getattr(exec_link, "login")
        try:
            login_success, expiry_datetime = func(cmd_pointer)

            if login_success == True:
                if not suppress:
                    print('You were successfully logged in.')
                return login_success, expiry_datetime
            else:
                if not suppress:
                    print('Login failed.')
                return False, None

            # write_login_registry(cmd_pointer.login_settings.copy())

        except BaseException as err:
            print(err)
            print("Login failed...")

            return False, None
