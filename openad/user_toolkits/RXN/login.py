"""Login Procedure for the RXN Toolkit"""

import os
import datetime
import importlib.util as ilu
from datetime import datetime, timezone
from rxn4chemistry import RXN4ChemistryWrapper
from openad.helpers.output import output_text, output_error, output_warning
from openad.helpers.output_msgs import msg
from openad.helpers.credentials import load_credentials, get_credentials, write_credentials
from openad.helpers.general import load_tk_module

API_CONFIG_BLANK = {"host": "None", "auth": {"username": "None", "api_key": "None"}, "verify_ssl": "false"}
DEFAULT_URL = "https://rxn.app.accelerate.science"

# Initialize the rxn client from the config file
# Input parameters for the example flow


def reset(cmd_pointer):
    """remove on reset signal, the underlying API key file to trigger reset"""
    cred_file = os.path.expanduser(f"{cmd_pointer.home_dir}/rxn_api.cred")
    if os.path.expanduser(cred_file):
        os.remove(cred_file)


def login(cmd_pointer):
    """logs onto the RXN service"""

    # Load module from toolkit folder
    rxn_helper = load_tk_module(cmd_pointer, "RXN", "rxn_include", "rxn_helper")()

    cred_file = os.path.expanduser(f"{cmd_pointer.home_dir}/rxn_api.cred")

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
                    return_val=False,
                )
            name, prj_id = rxn_helper.get_current_project(cmd_pointer)
        except Exception as e:  # pylint: disable=broad-exception-caught
            output_error(msg("err_login", "RXN", "Unable to connect to RXN Server"), return_val=False)
            output_error(msg("err_login", "RXN", f"system error {e}"), return_val=False)
            output_error(
                msg("if Persists,try removing RXN toolkit and restarting kernel or application./n Then add RXN again")
            )
            return False, None

        if name != cmd_pointer.settings["workspace"]:
            try:
                rxn_helper.sync_up_workspace_name(cmd_pointer)
            except:
                return False, None
        now = datetime.timestamp(now)
        return True, None  # No expiry on RXN Handles

    # if no Authentication file ask for authentication details and create
    try:
        config_file = get_creds(cred_file, cmd_pointer)
    except BaseException:
        return False, None
    x = cmd_pointer.login_settings["toolkits"].index("RXN")
    # fix for defaults when automating rxn cred application
    if config_file["host"].strip() == "None":
        config_file["host"] = ""

    try:
        client = RXN4ChemistryWrapper(api_key=config_file["auth"]["api_key"], base_url=config_file["host"])
        email = client.current_user()["response"]["payload"]["email"]
        if login_reset is True or first_login is True:
            workspace = cmd_pointer.settings["workspace"]
            output_text(
                f"<success>logging into RXN as: </success> {email}\n <success>Workspace: </success> {workspace}",
                return_val=False,
            )
        cmd_pointer.login_settings["toolkits_api"][x] = config_file["auth"]["api_key"]
        cmd_pointer.login_settings["client"][x] = client
        try:
            rxn_helper.sync_up_workspace_name(cmd_pointer)
        except:
            return False, None

        return True, None
    except Exception as e:  # pylint: disable=broad-exception-caught
        return False, None


def get_creds(cred_file, cmd_pointer):
    """get the nominated API key for the LLM"""
    api_config = load_credentials(cred_file)
    if api_config is None:
        output_warning("Please provide your RXN credentials:", return_val=False)
        output_text(
            f"<soft>Leave this blank to use the default: {DEFAULT_URL}</soft>",
            return_val=False,
        )
        api_config = API_CONFIG_BLANK.copy()
        api_config = get_credentials(
            cmd_pointer=cmd_pointer, credentials=api_config, creds_to_set=["host", "auth:api_key"]
        )
        write_credentials(api_config, cred_file)
    return api_config
