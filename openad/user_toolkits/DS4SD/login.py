"""Login Library for the Deepsearch Plugin"""
import os
import datetime
from datetime import datetime, timezone
import time
import requests
import jwt
import deepsearch as ds
from openad.helpers.output import msg, output_text, output_error, output_warning
from openad.helpers.credentials import load_credentials, get_credentials, write_credentials

DEFAULT_URL = "https://sds.app.accelerate.science/"
API_CONFIG_BLANK = {"host": "None", "auth": {"username": "None", "api_key": "None"}, "verify_ssl": "False"}


# Initialize the Deep Search client from the config file
# Input parameters for the example flow
def reset(cmd_pointer):
    """removes the deepsearch credentiuals file"""
    cred_file = os.path.expanduser(f"{cmd_pointer.home_dir}/deepsearch_api.cred")
    if os.path.isfile(cred_file):
        os.remove(cred_file)


def uri_valid(url: str) -> bool:
    """checks if a uri is valid"""
    try:
        request = requests.get(url, stream=True, timeout=10)
    except:  # pylint: disable=bare-except
        return False
    if request.status_code == 200:
        return True
    else:
        return False


def login(cmd_pointer):
    """Logs the Framework into the DS4SD service as defined by the user"""
    cred_file = os.path.expanduser(f"{cmd_pointer.home_dir}/deepsearch_api.cred")

    if not os.path.isfile(cred_file):
        login_reset = True
    else:
        login_reset = False

    first_login = False

    if "DS4SD" not in cmd_pointer.login_settings["toolkits"]:
        cmd_pointer.login_settings["toolkits"].append("DS4SD")
        cmd_pointer.login_settings["toolkits_details"].append({"type": "config_file", "session": "handle"})
        cmd_pointer.login_settings["toolkits_api"].append(None)
        cmd_pointer.login_settings["client"].append(None)
        cmd_pointer.login_settings["expiry"].append(None)
        x = cmd_pointer.login_settings["toolkits"].index("DS4SD")
        first_login = True
    elif login_reset is False:
        now = datetime.now(timezone.utc)
        x = cmd_pointer.login_settings["toolkits"].index("DS4SD")
        now = datetime.timestamp(now)
        expiry_time = cmd_pointer.login_settings["expiry"][x]

        if expiry_time is not None and expiry_time > now:
            expiry_datetime = time.strftime("%a %b %e, %G  at %R", time.localtime(expiry_time))
            return True, expiry_datetime

    cred_config = get_creds(cred_file, cmd_pointer)

    if cred_config["host"].strip() == "":
        cred_config["host"] = DEFAULT_URL
    if uri_valid(cred_config["host"]) is False:
        output_error(
            "Invalid URL Provided please check URL or VPN and try again", cmd_pointer=cmd_pointer, return_val=False
        )
        return False
    if cred_config["auth"]["username"].strip() == "":
        output_error("Invalid username provided try again", cmd_pointer=cmd_pointer, return_val=False)
        return False
    if cred_config["auth"]["api_key"].strip() == "":
        output_error("Invalid api key provided try again", cmd_pointer=cmd_pointer, return_val=False)
        return False

    try:
        x = cmd_pointer.login_settings["toolkits"].index("DS4SD")

        config = ds.DeepSearchConfig(host=cred_config["host"], verify_ssl=False, auth=cred_config["auth"])
        client = ds.CpsApiClient(config)
        api = ds.CpsApi(client)

        cmd_pointer.login_settings["toolkits_api"][x] = api
        cmd_pointer.login_settings["client"][x] = client

        cb = client.bearer_token_auth
        bearer = cb.bearer_token
        # decode jwt token
        decoded_token = jwt.decode(bearer, options={"verify_at_hash": False, "verify_signature": False}, verify=False)

        # Extract expiry time from token payload
        expiry_time = decoded_token["exp"]

        # Convert expiry time to a human-readable format
        expiry_datetime = time.strftime("%a %b %e, %G  at %R", time.localtime(expiry_time))
        if login_reset is True or first_login is True:
            workspace = cmd_pointer.settings["workspace"]
            email = cred_config["auth"]["username"]
            output_text(
                f"<success>loggining into DeepSearch as: </success> {email}\n <success>Workspace: </success> {workspace}",
                cmd_pointer=cmd_pointer,
                return_val=False,
            )
        cmd_pointer.login_settings["expiry"][x] = expiry_time

        return True, expiry_datetime
    except Exception as e:  # pylint: disable=broad-exception-caught
        output_error(
            msg("err_login", "DS4SD", f"Unable to connect to {config.host}", split=True),
            cmd_pointer=cmd_pointer,
            return_val=False,
        )
        return False, None


def get_creds(cred_file, cmd_pointer):
    """get the nominated API key for the LLM"""
    api_config = load_credentials(cred_file)
    if api_config is None:
        output_warning("Please Enter in Credentials for Deep Search", cmd_pointer=cmd_pointer, return_val=False)
        output_text(
            f"Enter the URL / Hostname: if the hostname is left blank it will default to '{DEFAULT_URL}' ",
            cmd_pointer=cmd_pointer,
            return_val=False,
        )
        api_config = API_CONFIG_BLANK.copy()
        api_config = get_credentials(
            cmd_pointer=cmd_pointer, credentials=api_config, creds_to_set=["host", "auth:username", "auth:api_key"]
        )
        write_credentials(api_config, cred_file)
    return api_config
