""" Provides Settings Management and Toolkit Management"""
import pickle
import os
import datetime
import tarfile
import shutil

# Core
from openad.core.grammar import create_statements

# Global variables
from openad.app.global_var_lib import _meta_dir_toolkits
from openad.app.global_var_lib import _meta_registry
from openad.app.global_var_lib import _meta_registry_session
from openad.app.global_var_lib import _meta_registry_settings

# Helpers
from openad.helpers.output import msg, output_text, output_error, output_success
from openad.helpers.general import confirm_prompt, refresh_prompt, other_sessions_exist


def initialise_registry():
    """initialise the registry/Settings File"""
    try:
        with open(_meta_registry, "wb") as handle:
            settings = pickle.dump(_meta_registry_settings, handle)
        return True
    except Exception as err:
        output_error("Error initialising the registry: " + err)
        return False


def update_main_registry_env_var(cmd_pointer, var_name, value):
    """update the registry file  wiht current version in memory"""
    registry = load_registry(cmd_pointer, orig_reg=True)
    registry["env_vars"][var_name] = value
    write_registry(registry, cmd_pointer, orig_reg=True)


# Load the user's registry data.


def load_registry(cmd_pointer, orig_reg=False):
    """load the registry file from disk
    if Orig_reg is false it loads the session settings file if not it saves
    back to the master registiry file"""
    if orig_reg is True:
        registry_file = _meta_registry
    else:
        registry_file = _meta_registry_session + cmd_pointer.session_id
    try:
        with open(registry_file, "rb") as handle:
            registry = pickle.loads(handle.read())
    except Exception:
        registry = _meta_registry_settings
    # New Section to test for updated Registry elements and add them
    # possibly is now redundant / Deprecated leaving for template
    # to show how we can manage settings file changes cross version
    try:
        dummy = registry["env_vars"]
    except Exception:
        registry["env_vars"] = {}

    return registry


# Write the user's registry data.
def write_registry(registry: dict, cmd_pointer, orig_reg=False):
    """Write the registry file from disk
    if Orig_reg is false it writes the session settings file if not it saves
    back to the master registiry file"""
    if orig_reg is True:
        registry_file = _meta_registry
    else:
        registry_file = _meta_registry_session + cmd_pointer.session_id
    try:
        with open(registry_file, "wb") as handle:
            settings = pickle.dump(registry, handle)
    except Exception as err:
        output_error("error writing registry" + err, cmd_pointer=cmd_pointer, return_val=False)
        return False
    return True


# Delete a individual session registry when the session ends.
def delete_session_registry(session_id: str):
    "deletes the session registry file based on session id"
    os.remove(_meta_registry_session + session_id)
    # os.remove('*' + session_id)
    return


# Determine when Other Sessions Exist


# Sledge Hammer remove other sessions, this enables new workspaces and toolkits to be created.
def clear_other_sessions(cmd_pointer, parser):
    """burned earth strategy for cleaning up Sessions directory"""
    if not confirm_prompt(msg("confirm_clear_sessions")):
        output_text(msg("no_sessions_cleared"), cmd_pointer, pad_btm=1, return_val=False)
        return

    file_list = os.listdir(os.path.dirname(_meta_registry_session))
    try:
        file_list.remove("registry.pkl" + cmd_pointer.session_id)
        for i in file_list:
            os.remove(os.path.dirname(_meta_registry_session) + "/" + i)
        # raise Exception('This is a test exception')
        return output_success(msg("success_clear_sessions"), cmd_pointer, pad_btm=1, return_val=False)
    except Exception as err:
        output_error(msg("err_clear_sessions", err, split=True), cmd_pointer, pad_btm=1, return_val=False)
        return False


def registry_add_toolkit(cmd_pointer, parser):
    """copies a toolkit from the installed source to the toolkits directory"""
    cmd_pointer.refresh_vector = True
    cmd_pointer.refresh_train = True
    cmd_pointer.settings["env_vars"]["refresh_help_ai"] = True
    update_main_registry_env_var(cmd_pointer, "refresh_help_ai", True)
    other, error_msg = other_sessions_exist(cmd_pointer)
    if other is True:
        return error_msg

    toolkit_name = parser["toolkit_name"]
    try:
        with open(_meta_registry, "rb") as handle:
            settings = pickle.loads(handle.read())
        if toolkit_name.upper() not in settings["toolkits"]:
            full_original_directory_name = (
                os.path.dirname(os.path.abspath(__file__)) + "/../user_toolkits/" + toolkit_name.upper() + "/"
            )
            target_directory = _meta_dir_toolkits + "/" + toolkit_name.upper()
            try:
                shutil.rmtree(target_directory + "/", ignore_errors=True)

            except Exception as err:
                output_error("Unable to write registry file::" + str(err), retun_val=False, cmd_pointer=False)

            shutil.copytree(full_original_directory_name, target_directory, dirs_exist_ok=True)
            settings["toolkits"].append(toolkit_name.upper())
            cmd_pointer.settings["toolkits"].append(toolkit_name.upper())
            write_registry(settings, cmd_pointer)  # moved this in scope was not being called hence
            write_registry(settings, cmd_pointer, True)

            if not os.path.isdir(_meta_dir_toolkits + "/" + toolkit_name.upper()):
                os.mkdir(_meta_dir_toolkits + "/" + toolkit_name.upper())
            else:
                datetime_str = str(datetime.datetime.now())
                archive_dir = f"{_meta_dir_toolkits}/ARCHIVE"
                archive_file = archive_dir + "/" + toolkit_name.upper() + "_" + datetime_str + ".tar.gz"

                if not os.path.isdir(archive_dir):
                    os.mkdir(archive_dir)

                with tarfile.open(archive_file, mode="w:gz") as archive:
                    archive.add(_meta_dir_toolkits + "/" + toolkit_name.upper(), recursive=True)
                    archive.close()

            # raise Exception('This is a test exception')
            return output_success(msg("success_toolkit_install", toolkit_name.upper()), cmd_pointer)
        else:
            return output_text(msg("toolkit_already_installed", toolkit_name.upper()), cmd_pointer, pad=1)

    except FileNotFoundError:
        return output_error(msg("invalid_toolkit", toolkit_name.upper(), split=True), cmd_pointer)

    except Exception as err:
        return output_error(msg("err_toolkit_install", toolkit_name.upper(), err, split=True), cmd_pointer)


def registry_deregister_toolkit(cmd_pointer, parser):
    """removes a toolkit from the registry"""
    other, error_msg = other_sessions_exist(cmd_pointer)
    cmd_pointer.refresh_vector = True
    cmd_pointer.refresh_train = True
    cmd_pointer.settings["env_vars"]["refresh_help_ai"] = True
    update_main_registry_env_var(cmd_pointer, "refresh_help_ai", True)
    if other is True:
        return error_msg

    toolkit_name = parser["toolkit_name"].upper()
    try:
        with open(_meta_registry, "rb") as handle:
            settings = pickle.loads(handle.read())
            if toolkit_name in settings["toolkits"]:
                settings["toolkits"].remove(toolkit_name.upper())
                cmd_pointer.settings["toolkits"].remove(toolkit_name.upper())
                write_registry(cmd_pointer.settings, cmd_pointer)

            else:
                return output_error(msg("fail_toolkit_not_registered", toolkit_name), cmd_pointer)

            if cmd_pointer.settings["context"] == toolkit_name:
                settings["context"] = None
                cmd_pointer.settings["context"] = None
                cmd_pointer.toolkit_current = None
                cmd_pointer.current_help.reset_help()
                write_registry(cmd_pointer.settings, cmd_pointer)
                refresh_prompt(cmd_pointer.settings)
                create_statements(cmd_pointer)

    except Exception as err:
        return output_error(msg("err_toolkit_remove", toolkit_name, err, split=True), cmd_pointer)
    write_registry(settings, cmd_pointer)
    write_registry(settings, cmd_pointer, True)
    return output_success(msg("success_toolkit_remove", toolkit_name.upper()), cmd_pointer)
