# Look for @Developer to find the places you need to edit.

import os
import jwt
import time
import requests
import datetime
from datetime import datetime, timezone
from openad.helpers.output import msg, output_text, output_error, output_warning, output_success
from openad.helpers.credentials import load_credentials, get_credentials, write_credentials


# @Developer Set the appropriate values:
TOOLKIT_NAME = "DEMO"
AUTH_URL_DEFAULT = "https://demo.com/auth"
PROMPT_AUTH_URL = False  # True = Prompt user for URL (default will be suggested.)

# Config template.
API_CONFIG_BLANK = {"host": AUTH_URL_DEFAULT, "auth": {"username": "None", "api_key": "None"}, "verify_ssl": "False"}


def login(cmd_pointer):
    """
    Authenticate the user with the toolkit's service.
    Triggered by `set context <toolkit_name>`
    """

    first_login = False

    # Try loading the credentials file.
    cred_file = os.path.expanduser(f"{cmd_pointer.home_dir}/{TOOLKIT_NAME.lower()}_api.cred")
    if not os.path.isfile(cred_file):
        login_reset = True
    else:
        login_reset = False

    # Prepare login settings when logging in for the first time.
    if TOOLKIT_NAME not in cmd_pointer.login_settings["toolkits"]:
        cmd_pointer.login_settings["toolkits"].append(TOOLKIT_NAME)
        cmd_pointer.login_settings["toolkits_details"].append({"type": "config_file", "session": "handle"})
        cmd_pointer.login_settings["toolkits_api"].append(None)
        cmd_pointer.login_settings["client"].append(None)
        cmd_pointer.login_settings["expiry"].append(None)
        first_login = True

    # Skip if you're already logged in during this session.
    elif login_reset is False:
        now = datetime.now(timezone.utc)
        now = datetime.timestamp(now)
        tk_index = cmd_pointer.login_settings["toolkits"].index(TOOLKIT_NAME)
        expiry_time = cmd_pointer.login_settings["expiry"][tk_index]

        if expiry_time is not None and expiry_time > now:
            return True, _readable_time(expiry_time)

    # Load credentials from file or prompt user.
    cred_config = _get_creds(cred_file, cmd_pointer)

    # Verify the auth URL is valid and responsive.
    valid, status = _uri_valid(cred_config["host"])
    if valid is False:
        if status == "invalid":
            output_error("Invalid host URL, please check for typos", cmd_pointer=cmd_pointer, return_val=False)
        else:
            output_error("Host URL unreachable, response code {status}", cmd_pointer=cmd_pointer, return_val=False)
        return False

    # Verify the username and API key are valid.
    if cred_config["auth"]["username"].strip() == "":
        output_error("Invalid username, please try again", cmd_pointer=cmd_pointer, return_val=False)
        return False
    if cred_config["auth"]["api_key"].strip() == "":
        output_error("Invalid API key, please try again", cmd_pointer=cmd_pointer, return_val=False)
        return False

    # Attempt login.
    try:
        # - - - @Developer start - - -
        # This section may need to be modified to fit your toolkit's API.
        client = None  # @Phil explain
        api = None  # @Phil explain
        auth_token = None  # <-- Output of this section
        # - - - @Developer end - - -

        # Get expiration time.
        expiry_time = jwt.decode(
            auth_token, options={"verify_at_hash": False, "verify_signature": False}, verify=False
        )["exp"]

        # Update login settings.
        tk_index = cmd_pointer.login_settings["toolkits"].index(TOOLKIT_NAME)
        cmd_pointer.login_settings["toolkits_api"][tk_index] = api
        cmd_pointer.login_settings["client"][tk_index] = client
        cmd_pointer.login_settings["expiry"][tk_index] = expiry_time

        # Success
        return True, _readable_time(expiry_time)

    # Fail
    except Exception as err:  # pylint: disable=broad-exception-caught
        return False, err


def reset(cmd_pointer):
    """
    Remove toolkit credentials file.
    Triggered by `set context <toolkit_name> reset`
    """
    cred_file = os.path.expanduser(f"{cmd_pointer.home_dir}/{TOOLKIT_NAME.lower()}_api.cred")
    if os.path.isfile(cred_file):
        os.remove(cred_file)


#
#
#


# Load credentials from pickle file, or prompt user is none present.
def _get_creds(cred_file, cmd_pointer):
    creds_config = load_credentials(cred_file)
    if creds_config is None:
        output_warning(f"Please provide your {TOOLKIT_NAME} credentials", cmd_pointer=cmd_pointer, return_val=False)
        if PROMPT_AUTH_URL:
            output_text(
                f"<soft>Leave this blank to use the default: {AUTH_URL_DEFAULT}</soft>",
                cmd_pointer=cmd_pointer,
                return_val=False,
            )
        creds_config = API_CONFIG_BLANK.copy()
        creds_config = get_credentials(
            cmd_pointer=cmd_pointer,
            credentials=creds_config,
            creds_to_set=["host", "auth:username", "auth:api_key"],
        )
        write_credentials(creds_config, cred_file)

    return creds_config


# Check if a uri is valid.
def _uri_valid(url: str) -> bool:
    try:
        request = requests.get(url, stream=True, timeout=10)
    except:  # pylint: disable=bare-except
        return False, "invalid"
    if request.status_code == 200:
        return True, request.status_code
    else:
        return False, request.status_code


# Convert timestamp to human-readable format.
def _readable_time(time_stamp):
    return time.strftime("%a %b %e, %G  at %R", time.localtime(time_stamp))
