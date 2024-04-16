import pyparsing as py
import os
import json
import glob
from openad.helpers.output import output_text, output_table, output_warning, output_error, output_success
from openad.openad_model_plugin.services import ModelService
from typing import List, Dict


SERVICE_DEFINTION_PATH = os.path.expanduser("~/.openad_model_services/")
SERVICES_PATH = os.path.expanduser("/definitions/services/")

model_service = ModelService()


def help_dict_create(
    name: str,  # Name of the comand - used for ...?
    command: str,  # Command structure, used for help, docs, training
    description: str,  # Description of the command, used for help, docs, training
    note: str = None,  # Additional note to the command, only used in help (eg. To learn more about runs, run `run ?`)
    url: str = None,  # Currently not used - URL to the documentation of the command?
    category: str = "Uncategorized",  # Category used to organize the commands in help & docs
    parent: str = None,  # Parent command, only relevant for follow-up commands like `result open`
):
    """Create a help dictionary"""
    return {
        "category": category,
        "name": name,
        "command": command,
        "description": description,
        "note": note,
        "url": url,
        "parent": parent,
    }


def get_service_defs(reference) -> list:
    """pulls the list of available service definitions"""
    service_list = []
    service_files = glob.glob(reference + "/*.json")
    for file in service_files:
        with open(file, "r") as file_handle:
            try:
                jdoc = json.load(file_handle)
                service_list.append(jdoc)
            except Exception as e:
                print(e)
                print("invalid service json definition  " + file)
    return service_list


def get_cataloged_service_defs():
    """Returns a list of cataloged Services definitions and their Namespaces"""
    if not os.path.exists(SERVICE_DEFINTION_PATH):
        os.makedirs(SERVICE_DEFINTION_PATH)

    list_of_namespaces = [
        os.path.basename(f.path) for f in os.scandir(SERVICE_DEFINTION_PATH) if f.is_dir()
    ]  # os.walk(SERVICE_DEFINTION_PATH)

    service_list_by_catalog = {}
    for namespace in list_of_namespaces:
        service_list = []
        services_path = SERVICE_DEFINTION_PATH + namespace + SERVICES_PATH
        if os.path.exists(services_path):
            service_list = get_service_defs(services_path)
        service_list_by_catalog[namespace] = service_list
    return service_list_by_catalog


def service_status(service_name) -> Dict:
    """Get a model service status"""
    return json.loads(model_service.status(service_name))


def service_up(service_name) -> None:
    """This function synchronously starts a service"""
    model_service.up(service_name)


def service_down(service_name) -> None:
    """This function synchronously shuts down a service"""
    model_service.down(service_name)


def remove_cataloged_service(service_name) -> None:
    """This function removes a catalog from the ~/.openad_model_service directory"""
    model_service.remove_service(service_name)


def get_service_endpoint(service_name) -> List:
    """gets the service endpoint for a given service, if endpoint is not available it returns None"""
    if service_name is None:
        "may in future return a default local service"
        return []
    endpoint = json.loads(model_service.status(service_name)).get("url")
    return endpoint


def service_catalog_grammar(statements: list, help: list, service_list: list):
    """This function creates the required grammar for managing cataloging services and model up or down"""
    add = py.CaselessKeyword("add")
    model = py.CaselessKeyword("model")
    service = py.CaselessKeyword("service")
    fr_om = py.CaselessKeyword("from")
    path = py.CaselessKeyword("path")
    quoted_string = py.QuotedString("'", escQuote="\\")
    a_s = py.CaselessKeyword("as")

    statements.append(
        py.Forward(add + model + service + fr_om + path + quoted_string("path") + a_s + quoted_string("service_name"))(
            "add_model_path"
        )
    )
    help.append(
        help_dict_create(
            name="Add Model from Path",
            category="General",
            command="add model from path '<path to model directory>' as '<service_name>'",
            description="add a model definition to the catalog.",
        )
    )
