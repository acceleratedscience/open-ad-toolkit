import os
import pickle
from typing import Dict, TypedDict
from openad.openad_model_plugin.utils import get_logger, bcolors
from openad.openad_model_plugin.config import SERVICE_MODEL_PATH


logger = get_logger(__name__, color=bcolors.OKCYAN + bcolors.UNDERLINE)


# path to authentication lookup table
auth_lookup_path = os.path.join(SERVICE_MODEL_PATH, "auth_lookup.pkl")


class LookupTable(TypedDict):
    """helper class for python types"""

    auth_table: Dict[str, str]
    service_table: Dict[str, str]


def load_lookup_table() -> LookupTable:
    """load lookup authentication lookup table from pickle file"""
    logger.debug("loading auth lookup table")
    # load the current auth groups
    if os.path.exists(auth_lookup_path):
        with open(auth_lookup_path, "rb") as file:
            data = pickle.load(file)
    else:
        data = {"auth_table": {}, "service_table": {}}
        with open(auth_lookup_path, "wb") as file:
            pickle.dump(data, file)
    return data


def update_lookup_table(group_name, api_key=None, service=None):
    """update the lookup table values on either api_key or model service"""
    logger.debug(f"updating auth group '{group_name}' {api_key=} {service=}")
    # Load the existing data or create a new dictionary
    data = load_lookup_table()
    # Update the dictionary with new key-value pairs
    if api_key:
        # create auth group entry
        data["auth_table"].update({group_name: api_key})
    if service and group_name in data["auth_table"]:
        # map a model service to auth group
        data["service_table"].update({service: group_name})
    # Save the updated dictionary back to the pickle file
    with open(auth_lookup_path, "wb") as file:
        pickle.dump(data, file)
    # return latest data
    return data


def remove_service_group(service_name):
    """remove a model service from lookup table"""
    logger.debug(f"removing service auth group | '{service_name}'")
    # Load the existing data or create a new dictionary
    data = load_lookup_table()
    data["service_table"].pop(service_name)
    # Save the updated dictionary back to the pickle file
    with open(auth_lookup_path, "wb") as file:
        pickle.dump(data, file)
    # return latest data
    return data


def remove_auth_group(group_name):
    """remove a authentication group from lookup table and and relationships to model service"""
    logger.debug(f"removing auth group '{group_name}'")
    # Load the existing data or create a new dictionary
    data = load_lookup_table()
    # remove the dictionary values with matching auth group
    data["auth_table"].pop(group_name)
    for name, val in data["service_table"].copy().items():
        logger.debug(f"evaluating {name=} {val=}")
        if group_name == val:
            logger.debug(f"removing auth group '{group_name}' from service '{name}'")
            data["service_table"].pop(name)
    # Save the updated dictionary back to the pickle file
    with open(auth_lookup_path, "wb") as file:
        pickle.dump(data, file)
    # return latest data
    return data
