import pyparsing as py
import os
import json
import glob
from openad.helpers.output import output_text, output_table, output_warning, output_error, output_success

SERVICE_DEFINTION_PATH = os.path.expanduser("~/.openad_model_services/")
SERVICES_PATH = os.path.expanduser("/definitions/services/")


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


def get_services(reference) -> list:
    """pulls the list of available services for"""

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


def get_cataloged_services():
    """Returns a list of cataloged Services and their Namespaces"""
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
            service_list = get_services(SERVICE_DEFINTION_PATH + namespace + SERVICES_PATH)

        service_list_by_catalog[namespace] = service_list

    return service_list_by_catalog


def list_cataloged_services(cmd_pointer, parser):
    """This function catalogs a service"""
    pass


def catalog_service(cmd_pointer, parser):

    pass


def service_up(cmd_pointer, parser):
    """This function synchronously starts a service"""
    pass


def service_down(cmd_pointer, parser):
    """This function synchronously shuts down a service"""
    pass


def remove_cataloged_service(cmd_pointer, parser):
    """This function removes a catalog from the ~/.openad_model_service directory"""
    pass


def get_service_endpoint(service_name):
    """gets the service endpoint for a given service, if endpoint is not available it returns None"""
    Endpoint = None
    if service_name is None:
        "may in future return a default local service"
        return Endpoint
    # this is a mock-up while service API is integrated
    if service_name == "gt4sd_gen":
        Endpoint = "http://127.0.0.1:8090"
    elif service_name == "gt4sd_prop":
        Endpoint = "http://127.0.0.1:8080"

    return Endpoint


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
