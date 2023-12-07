""" Helpers for getting and Saving Credentials"""
import pickle
import readline
import os
from openad.helpers.general import user_input, user_secret

DEFAULT_SECRET = ["api_key", "password"]
DEFAULT_USER = ["user_name", "email", "name"]
DEFAULT_CREDENTIALS = {"host": "None", "auth": {"user_name": "None", "api_key": "None"}, "verify_ssl": "false"}
DEFAULT_CREDS_TO_SET = ["host", "auth:user_name", "auth:api_key"]


def get_credentials(cmd_pointer, credentials=DEFAULT_CREDENTIALS, creds_to_set=DEFAULT_CREDS_TO_SET) -> dict:
    """prompts user for set of credentials"""
    new_credentials = credentials.copy()
    for cred in creds_to_set:
        value_to_prompt = cred.split(":", maxsplit=1)[-1]
        if value_to_prompt in DEFAULT_SECRET:
            prompted_value = user_secret(cmd_pointer, value_to_prompt)
        elif value_to_prompt in DEFAULT_USER:
            prompted_value = user_input(cmd_pointer, value_to_prompt)
        else:
            prompted_value = user_input(cmd_pointer, value_to_prompt)
        if readline.get_current_history_length() > 0:
            readline.remove_history_item(readline.get_current_history_length() - 1)
        new_credentials = assign_dict_value(new_credentials, cred, prompted_value)
    return new_credentials


def assign_dict_value(input_dict: dict, path: str, value) -> dict:
    "assigns a value to a defined path in a dictionary"
    output_dict = input_dict.copy()
    statement = "exec('output_dict"
    for key in path.split(":"):
        statement = f'{statement}["{key}"]'
    if isinstance(value, (str)):
        statement = statement + f'="{value}"'
    else:
        statement = statement + f"={value}"
    eval(statement + "')", {"output_dict": output_dict})  # pylint:  disable=eval-used
    return output_dict


def load_credentials(location):
    """Load the user's api login data."""
    if not os.path.exists(os.path.expanduser(location)):
        return None
    with open(os.path.expanduser(location), "rb") as handle:
        credentials = pickle.loads(handle.read())
        handle.close()
        return credentials


def write_credentials(registry: dict, location):
    """Dumps the LLM api details to home directory of app"""
    with open(os.path.expanduser(location), "wb") as handle:
        pickle.dump(registry, handle)
    return True


if __name__ == "__main__":
    print(load_credentials("~/.openad/openai_api.json"))
