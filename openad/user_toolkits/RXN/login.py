"""Login Procedure for the RXN Toolkit"""

import os
import datetime
import importlib.util as ilu
from datetime import datetime, timezone
from rxn4chemistry import RXN4ChemistryWrapper
from openad.helpers.output import msg, output_text, output_error, output_warning
from openad.helpers.credentials import load_credentials, get_credentials, write_credentials

API_CONFIG_BLANK = {"host": "None", "auth": {"username": "None", "api_key": "None"}, "verify_ssl": "false"}
DEFAULT_URL = "https://rxn.app.accelerate.science"

# Initialize the rxn client from the config file
# Input parameters for the example flow


def get_include_lib(cmd_pointer):
    """load rxn include library"""
    folder = cmd_pointer.toolkit_dir + "/RXN/rxn_include.py"
    file = "rxn_include"
    spec = ilu.spec_from_file_location(file, folder)
    rxn = ilu.module_from_spec(spec)
    spec.loader.exec_module(rxn)
    rxn_helper = rxn.rxn_helper()
    return rxn_helper


def reset(cmd_pointer):
    """remove on reset signal, the underlying API key file to trigger reset"""
    cred_file = os.path.expanduser(f"{cmd_pointer.home_dir}/rxn_api.cred")
    if os.path.expanduser(cred_file):
        os.remove(cred_file)


def login(cmd_pointer):
    """logs onto the RXN service"""
    cred_file = os.path.expanduser(f"{cmd_pointer.home_dir}/rxn_api.cred")
    rxn_helper = get_include_lib(cmd_pointer)

    if not os.path.isfile(cred_file):
        login_reset = True
    else:
        login_reset = False

    first_login = False  # Used primarily for Notebook mode to signal logged on only once

    # If we have not already logged onto the RXN service in this Session
    if "RXN" not in cmd_pointer.login_settings["toolkits"]:
        first_login = True
        cmd_pointer.login_settings["toolkits"].append("RXN")
        cmd_pointer.login_settings["toolkits_details"].append({"type": "config_file", "session": "handle"})
        cmd_pointer.login_settings["toolkits_api"].append(None)
        cmd_pointer.login_settings["client"].append(None)
        cmd_pointer.login_settings["expiry"].append(None)
        x = cmd_pointer.login_settings["toolkits"].index("RXN")
        cmd_pointer.login_settings["session_vars"].append(
            rxn_helper._RXN_VARS_TEMPLATE
        )  # pylint: disable=protected-access
    elif login_reset is False:  # If a login reset has been issued or there is no authentication file
        now = datetime.now(timezone.utc)
        x = cmd_pointer.login_settings["toolkits"].index("RXN")
        client = cmd_pointer.login_settings["client"][cmd_pointer.login_settings["toolkits"].index("RXN")]
        try:
            if login_reset is True or first_login is True:
                email = client.current_user()["response"]["payload"]["email"]
                workspace = cmd_pointer.settings["workspace"]
                output_text(
                    f"<success>logging into RXN as: </success> {email}\n <success>Workspace: </success> {workspace}",
                    cmd_pointer=cmd_pointer,
                    return_val=False,
                )
            name, prj_id = rxn_helper.get_current_project(cmd_pointer)
        except Exception as e:  # pylint: disable=broad-exception-caught
            output_error(
                msg("err_login", "RXN", "Unable to connect to RXN Server", split=True),
                cmd_pointer=cmd_pointer,
                return_val=False,
            )
            output_error(
                msg("err_login", "RXN", f"system error {e}", split=True), cmd_pointer=cmd_pointer, return_val=False
            )
            return False, None

        if name != cmd_pointer.settings["workspace"]:
            rxn_helper.sync_up_workspace_name(cmd_pointer)
        now = datetime.timestamp(now)
        return True, None  # No expiry on RXN Handles

    # if no Authentication file ask for authentication details and create
    config_file = get_creds(cred_file, cmd_pointer)

    x = cmd_pointer.login_settings["toolkits"].index("RXN")
    try:
        client = RXN4ChemistryWrapper(api_key=config_file["auth"]["api_key"], base_url=config_file["host"])
        email = client.current_user()["response"]["payload"]["email"]
        if login_reset is True or first_login is True:
            workspace = cmd_pointer.settings["workspace"]
            output_text(
                f"<success>logging into RXN as: </success> {email}\n <success>Workspace: </success> {workspace}",
                cmd_pointer=cmd_pointer,
                return_val=False,
            )
        cmd_pointer.login_settings["toolkits_api"][x] = config_file["auth"]["api_key"]
        cmd_pointer.login_settings["client"][x] = client
        rxn_helper.sync_up_workspace_name(cmd_pointer, reset=True)
        return True, None
    except Exception as e:  # pylint: disable=broad-exception-caught
        output_error(
            msg("err_login", "RXN", f"Unable to connect to  {config_file['host']}", split=True),
            cmd_pointer=cmd_pointer,
            return_val=False,
        )
        output_error(
            msg("err_login", "RXN", f"system error {e}", split=True), cmd_pointer=cmd_pointer, return_val=False
        )
        return False, None


def get_creds(cred_file, cmd_pointer):
    """get the nominated API key for the LLM"""
    api_config = load_credentials(cred_file)
    if api_config is None:
        output_warning("Setting Authentication Details for RXN:", cmd_pointer=cmd_pointer, return_val=False)
        output_text(
            f"Enter the Hostname: if the hostname is left blank it will default to '{DEFAULT_URL}' ",
            cmd_pointer=cmd_pointer,
            return_val=False,
        )
        api_config = API_CONFIG_BLANK.copy()
        api_config = get_credentials(
            cmd_pointer=cmd_pointer, credentials=api_config, creds_to_set=["host", "auth:api_key"]
        )
        write_credentials(api_config, cred_file)
    return api_config
