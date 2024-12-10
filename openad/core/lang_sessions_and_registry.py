""" Provides Settings Management and Toolkit Management"""

import pickle
import os
import datetime
import tarfile
import shutil
from IPython.display import display

# Core
from openad.core.grammar import create_statements
from openad.openad_model_plugin.utils import bcolors, get_logger

# Global variables
from openad.app.global_var_lib import _meta_dir_toolkits
from openad.app.global_var_lib import _meta_registry
from openad.app.global_var_lib import _meta_registry_session
from openad.app.global_var_lib import _meta_registry_settings
from openad.app.global_var_lib import _all_toolkits

# Helpers
from openad.helpers.output import output_text, output_error, output_success, output_warning
from openad.helpers.output_msgs import msg
from openad.helpers.general import confirm_prompt, refresh_prompt, other_sessions_exist

logger = get_logger(__name__, color=bcolors.OKCYAN + bcolors.UNDERLINE)


def initialise_registry():
    """
    Initialise the registry file.
    The registry file is a pickle file that contains the user's settings.
    """
    try:
        with open(_meta_registry, "wb") as handle:
            settings = pickle.dump(_meta_registry_settings, handle)
        return True
    except Exception as err:
        output_error("Error initialising the registry: " + err)
        return False


def update_main_registry_env_var(cmd_pointer, var_name, value):
    """
    Update the registry file with the current version in memory.
    """
    registry = load_registry(cmd_pointer, orig_reg=True)
    registry["env_vars"][var_name] = value
    write_registry(registry, cmd_pointer, orig_reg=True)


# Load the user's registry data.


def load_registry(cmd_pointer, orig_reg=False):
    """
    Load the registry file from disk.
    If orig_reg=False, it loads the session registry file,
    otherwise it loads the master registry file.
    """
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


def write_registry(registry: dict, cmd_pointer, orig_reg=False):
    """
    Write the registry file to disk.
    If orig_reg=False it writes to the session registry file,
    otherwise it writes to the master registry file.
    """
    if len(os.listdir(os.path.dirname(_meta_registry_session))) > 1 and orig_reg == True:
        logger.warning(
            "Main Registry was not Saved as other sessions are open. use clear Session command to clear zombie sessions"
        )
        return True

    if orig_reg is True:
        registry_file = _meta_registry
    else:
        registry_file = _meta_registry_session + cmd_pointer.session_id
    try:
        with open(registry_file, "wb") as handle:
            settings = pickle.dump(registry, handle)
    except Exception as err:
        output_error("error writing registry" + err, return_val=False)
        return False
    return True


def delete_session_registry(session_id: str):
    """
    Deletes the session registry file when the session ends.
    """
    try:
        os.remove(_meta_registry_session + session_id)
    except Exception as err:  # pylint: disable=broad-except
        # When clear sessions was called from another session,
        # the session registry file won't exist and we don't
        # want to throw an error.
        pass

    # os.remove('*' + session_id)
    return


def clear_sessions(cmd_pointer, parser):
    """
    Burned-earth approach for cleaning up the session directory.
    Required to enable new workspaces and toolkits to be created.
    """
    if not confirm_prompt(msg("confirm_clear_sessions")):
        output_error(msg("no_sessions_cleared"), return_val=False)
        return False

    file_list = os.listdir(os.path.dirname(_meta_registry_session))
    try:
        file_list.remove("registry.pkl" + cmd_pointer.session_id)
        for i in file_list:
            os.remove(os.path.dirname(_meta_registry_session) + "/" + i)
        # raise Exception("This is a test exception")
        return output_success(msg("success_clear_sessions"), return_val=False)
    except Exception as err:  # pylint: disable=broad-except
        output_error(msg("err_clear_sessions", err), return_val=False)
        return False


def registry_add_toolkit(cmd_pointer, parser, switch_context=True, suppress_output=False):
    """
    Installs a toolkit by copying it from the installation
    source (openad/user_toolkits) to the toolkits directory.
    """
    from openad.app.main_lib import set_context_by_name

    other_sesh = other_sessions_exist(cmd_pointer)
    if other_sesh is True:
        return False

    cmd_pointer.refresh_vector = True
    cmd_pointer.refresh_train = True
    cmd_pointer.settings["env_vars"]["refresh_help_ai"] = True
    update_main_registry_env_var(cmd_pointer, "refresh_help_ai", True)

    toolkit_name = parser["toolkit_name"]
    try:
        # raise Exception("This is a test exception")
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
                # if not suppress_output:
                output_error("Unable to write registry file::" + str(err), retun_val=False)
                return False

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

            # Set context to new toolkit
            if switch_context:
                set_context_by_name(cmd_pointer, toolkit_name.upper())

            if not suppress_output:
                output_success(msg("success_toolkit_install", toolkit_name.upper()), return_val=False)
            return True
        else:
            if not suppress_output:
                output_warning(msg("toolkit_already_installed", toolkit_name.upper()), return_val=False)
            return False

    except FileNotFoundError:
        # if not suppress_output:
        output_error(msg("invalid_toolkit", toolkit_name.upper()), return_val=False)
        return False

    except Exception as err:
        # if not suppress_output:
        output_error(msg("err_toolkit_install", toolkit_name.upper(), err), return_val=False)
        return False


def registry_remove_toolkit(cmd_pointer, parser, suppress_output=False):
    """
    Removes a toolkit from the registry.
    """
    other_sesh = other_sessions_exist(cmd_pointer)
    if other_sesh is True:
        return False

    cmd_pointer.refresh_vector = True
    cmd_pointer.refresh_train = True
    cmd_pointer.settings["env_vars"]["refresh_help_ai"] = True
    update_main_registry_env_var(cmd_pointer, "refresh_help_ai", True)

    toolkit_name = parser["toolkit_name"].upper()
    try:
        # raise Exception("This is a test exception")
        with open(_meta_registry, "rb") as handle:
            settings = pickle.loads(handle.read())
            if toolkit_name in settings["toolkits"]:
                settings["toolkits"].remove(toolkit_name.upper())
                cmd_pointer.settings["toolkits"].remove(toolkit_name.upper())
                write_registry(cmd_pointer.settings, cmd_pointer)

            else:
                # if not suppress_output:
                output_error(msg("fail_toolkit_not_registered", toolkit_name), return_val=False)
                return False

            if cmd_pointer.settings["context"] == toolkit_name:
                settings["context"] = None
                cmd_pointer.settings["context"] = None
                cmd_pointer.toolkit_current = None
                cmd_pointer.current_help.reset_help()
                write_registry(cmd_pointer.settings, cmd_pointer)
                refresh_prompt(cmd_pointer.settings)
                create_statements(cmd_pointer)

    except Exception as err:
        # if not suppress_output:
        output_error(msg("err_toolkit_remove", toolkit_name, err), return_val=False)
        return False

    write_registry(settings, cmd_pointer)
    write_registry(settings, cmd_pointer, True)
    if not suppress_output:
        output_success(msg("success_toolkit_remove", toolkit_name.upper()), return_val=False)
    return True


def update_toolkit(cmd_pointer, parser, suppress_output=False):
    """
    Updates a toolkit by removing and re-installing it.
    """
    from openad.app.main_lib import set_context_by_name

    # Get current context, so we can reset it in case
    # the current toolkit is updated.
    prev_toolkit = cmd_pointer.settings["context"]

    # Check if toolkit is installed
    toolkit_name = parser["toolkit_name"].upper()
    is_installed = bool(toolkit_name in cmd_pointer.settings["toolkits"])

    # Remove
    if is_installed:
        remove_success = registry_remove_toolkit(cmd_pointer, parser, suppress_output=True)
    else:
        remove_success = True

    # Re-install if removal was successful
    if remove_success:
        install_success = registry_add_toolkit(cmd_pointer, parser, switch_context=False, suppress_output=True)

        if install_success:
            # If the current toolkit was updated, we need
            # to reset the context to the current toolkit.
            if cmd_pointer.settings["context"] is not prev_toolkit:
                set_context_by_name(cmd_pointer, prev_toolkit, suppress_splash=True)

            # Display success message
            if not suppress_output:
                output_success(msg("success_update_toolkit", parser["toolkit_name"].upper()), return_val=False)
            return True


def update_all_toolkits(cmd_pointer, parser):
    """
    Updates all installed toolkits at once.
    """

    other_sesh = other_sessions_exist(cmd_pointer)
    if other_sesh is True:
        return

    update_successes = []
    update_fails = []

    # Assemble list of installed toolkits.
    installed_toolkits = []
    for toolkit_name in _all_toolkits:
        if toolkit_name in cmd_pointer.settings["toolkits"]:
            installed_toolkits.append(toolkit_name)

    # Update each toolkit.
    for toolkit_name in installed_toolkits:
        parser["toolkit_name"] = toolkit_name
        success = update_toolkit(cmd_pointer, parser, suppress_output=True)
        if success:
            update_successes.append(toolkit_name)
        else:
            update_fails.append(toolkit_name)

    # Display successes.
    if len(update_successes):
        list_success = []
        for toolkit_name in update_successes:
            list_success.append(f"<reset>- </reset><yellow>{toolkit_name}</yellow>")
        output_success(msg("success_update_all_toolkits", "\n".join(list_success)), return_val=False)

    # Display failures.
    if len(update_fails):
        list_fail = []
        for toolkit_name in update_fails:
            list_fail.append(f"<reset>- </reset><yellow>{toolkit_name}</yellow>")
        output_error(msg("err_update_all_toolkits", "\n".join(list_fail)), return_val=False)
