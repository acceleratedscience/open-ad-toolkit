import os
import pickle
from typing import Dict, TypedDict

from filelock import FileLock

from openad.openad_model_plugin.config import SERVICE_MODEL_PATH
from openad.openad_model_plugin.utils import bcolors, get_logger

logger = get_logger(__name__, color=bcolors.OKGREEN)


# path to authentication lookup table
auth_lookup_path = os.path.join(SERVICE_MODEL_PATH, "auth_lookup.pkl")
# lock file to prevent concurrent access
auth_lookup_lock = FileLock(os.path.join(SERVICE_MODEL_PATH, "auth_lookup.pkl.lock"))


class LookupTable(TypedDict):
    """helper class for python types"""

    auth_table: Dict[str, str]
    service_table: Dict[str, str]


def save_lookup_table(data: LookupTable):
    """save authentication lookup table to pickle file"""
    logger.debug("saving auth lookup table")
    with auth_lookup_lock:  # lock file to prevent concurrent access
        with open(auth_lookup_path, "wb") as file:
            pickle.dump(data, file)


def hide_api_keys(data: LookupTable) -> LookupTable:
    """hide api key from output. in place operation"""
    logger.debug("hiding api keys")
    for auth_group, api_key in data["auth_table"].items():
        data["auth_table"][auth_group] = api_key[:6] + "..." if len(api_key) > 4 else api_key
    return data


def load_lookup_table(hide_api: bool = False) -> LookupTable:
    """load authentication lookup table from pickle file"""
    logger.debug("loading auth lookup table")
    with auth_lookup_lock:  # lock file to prevent concurrent access
        # load the current auth groups
        if os.path.exists(auth_lookup_path):
            with open(auth_lookup_path, "rb") as file:
                data = pickle.load(file)
        else:
            data = {"auth_table": {}, "service_table": {}}
            save_lookup_table(data)  # create an empty lookup table
    # make only 6 character of api key visible
    if hide_api:
        logger.debug("hiding api keys1")
        return hide_api_keys(data)
    return data


def get_service_api_key(service_name: str) -> str:
    """get api key from auth lookup table. returns empty string for no api key"""
    # get lookup table
    auth_lookup_table = load_lookup_table()
    # find group name belonging to service
    auth_group = auth_lookup_table["service_table"].get(service_name, "")
    api_key = auth_lookup_table["auth_table"].get(auth_group, "")
    logger.debug(f"get service api key | {service_name=} {auth_group=} {api_key=}")
    return api_key


def update_lookup_table(auth_group, api_key=None, service=None, hide_api=False) -> LookupTable:
    """update the lookup table values on either api_key or model service"""
    logger.debug(f"updating auth group '{auth_group}' {api_key=} {service=}")
    # Load the existing data or create a new dictionary
    data = load_lookup_table()
    # Update the dictionary with new key-value pairs
    if api_key:
        # create auth group entry
        data["auth_table"].update({auth_group: api_key})
    if service and auth_group in data["auth_table"]:
        # map a model service to auth group
        data["service_table"].update({service: auth_group})
    # Save the updated dictionary back to the pickle file
    save_lookup_table(data)
    # return latest data
    if hide_api:
        return hide_api_keys(data)
    return data


def remove_service_group(service_name, hide_api=False) -> LookupTable:
    """remove a model service from lookup table"""
    logger.debug(f"removing service '{service_name}' from auth group lookup table")
    # Load the existing data or create a new dictionary
    data = load_lookup_table()
    data["service_table"].pop(service_name, None)
    # Save the updated dictionary back to the pickle file
    save_lookup_table(data)
    # return latest data
    if hide_api:
        return hide_api_keys(data)
    return data


def remove_auth_group(auth_group, hide_api=False) -> LookupTable:
    """remove a authentication group from lookup table and and relationships to model service"""
    logger.debug(f"removing auth group '{auth_group}'")
    # Load the existing data or create a new dictionary
    data = load_lookup_table()
    # remove the dictionary values with matching auth group
    data["auth_table"].pop(auth_group, None)
    for name, val in data["service_table"].copy().items():
        logger.debug(f"evaluating {name=} {val=}")
        if auth_group == val:
            logger.debug(f"removing auth group '{auth_group}' from service '{name}'")
            data["service_table"].pop(name, None)
    # Save the updated dictionary back to the pickle file
    save_lookup_table(data)
    # return latest data
    if hide_api:
        return hide_api_keys(data)
    return data
