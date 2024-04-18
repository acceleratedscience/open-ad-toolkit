import pyparsing as py
import os
import json
import glob
from openad.helpers.output import output_text, output_table, output_warning, output_error, output_success
from openad.openad_model_plugin.services import ServiceManager, UserProvidedConfig
from typing import List, Dict


SERVICE_DEFINTION_PATH = os.path.expanduser("~/.openad_model_services/")
SERVICES_PATH = os.path.expanduser("/definitions/services/")
if not os.path.exists(SERVICE_DEFINTION_PATH):
    os.makedirs(SERVICE_DEFINTION_PATH)


Service = ServiceManager()
### example of how to load services by namespace ###
# with Service(TEST_PATH) as service:
#     print(service.list())

def get_namespaces():
    list_of_namespaces = [
        os.path.basename(f.path) for f in os.scandir(SERVICE_DEFINTION_PATH) if f.is_dir()
    ]  # os.walk(SERVICE_DEFINTION_PATH)
    return list_of_namespaces

# load available model services on startup
try:
    for namespace in get_namespaces():
        with Service(SERVICE_DEFINTION_PATH + namespace) as service:
            print(f"> loaded: {service.cache}")
except:
    pass

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


def get_catalog_namespaces(cmd_pointer, parser) -> Dict:
    """Get a model service status"""
    ns = get_namespaces()
    print(ns)

def model_service_status(cmd_pointer, parser):
    """This function catalogs a service"""
    # get list of directory names for the catalog models
    namespaces = get_namespaces()
    ns_status = []
    for ns in namespaces:
        try:
            with Service(SERVICE_DEFINTION_PATH + ns) as model:
                    ns_status.append({ns: model.get_short_status(ns)}) # check list first
        except Exception as e:
            print(e)
            continue  # service not cataloged or doesnt exist
    print(ns_status)
    # print(json.dumps(ns_status[0], indent=2))

def catalog_add_model_service(cmd_pointer, parser):
    """Add model service repo to catalog"""
    service_name = parser.as_dict()["service_name"]
    path = parser.as_dict()["path"]
    # add service to api
    model_path = SERVICE_DEFINTION_PATH + service_name
    with Service(model_path) as model:
        # configure the sky yaml
        config = UserProvidedConfig(
            workdir=model_path,
            port=8080,
            setup="docker buildx build -f Dockerfile -t inference-service .",
            run="docker run --rm \
                    -p 8080:8080 \
                    inference-service",
            disk_size=100,
            )
        model.add_service(service_name, config)
        print(f"service {service_name} added to catalog")


def service_up(cmd_pointer, parser) -> None:
    """This function synchronously starts a service"""
    service_name = parser.as_dict()["service_name"]
    with Service(SERVICE_DEFINTION_PATH + service_name) as model:
        model.up(service_name)
        print(f"{service_name} started")


def service_down(cmd_pointer, parser) -> None:
    """This function synchronously shuts down a service"""
    service_name = parser.as_dict()["service_name"]
    with Service(SERVICE_DEFINTION_PATH + service_name) as model:
        model.down(service_name)
        print(f"{service_name} terminating")


def uncatalog_model_service(cmd_pointer, parser):
    """This function removes a catalog from the ~/.openad_model_service directory"""
    service_name = parser.as_dict()["service_name"]
    with Service(SERVICE_DEFINTION_PATH + service_name) as model:
        try:
            model.down(service_name)
        except Exception as e:
            pass
        finally:
            model.remove_service(service_name)
    print(f"service {service_name} removed from catalog")


def get_service_endpoint(service_name) -> List:
    """gets the service endpoint for a given service, if endpoint is not available it returns None"""
    if service_name is None:
        "may in future return a default local service"
        return []
    with Service(SERVICE_DEFINTION_PATH + service_name) as model:
        endpoint = json.loads(model.status(service_name)).get("url")
        return endpoint



def service_catalog_grammar(statements: list, help: list):
    """This function creates the required grammar for managing cataloging services and model up or down"""
    catalog = py.CaselessKeyword("catalog")
    uncatalog = py.CaselessKeyword("uncatalog")
    model = py.CaselessKeyword("model")
    up = py.CaselessKeyword("up")
    down = py.CaselessKeyword("down")
    service = py.CaselessKeyword("service")
    status = py.CaselessKeyword("status")
    fr_om = py.CaselessKeyword("from")
    _list = py.CaselessKeyword("list")
    path = py.CaselessKeyword("path")
    quoted_string = py.QuotedString("'", escQuote="\\")
    a_s = py.CaselessKeyword("as")

    statements.append(py.Forward(model + service + status)("model_service_status"))
    help.append(
        help_dict_create(
            name="model service status",
            category="Model",
            command="model service status",
            description="get the status of currently cataloged services",
        )
    )

    statements.append(py.Forward(model + catalog + _list)("get_catalog_namespaces"))
    help.append(
        help_dict_create(
            name="model catalog list",
            category="Model",
            command="model catalog list",
            description="get the list of currently cataloged services",
        )
    )

    statements.append(
        py.Forward(uncatalog + model + service + quoted_string("service_name"))("uncatalog_model_service")
    )
    help.append(
        help_dict_create(
            name="uncatalog model service",
            category="Model",
            command="uncatalog model service '<service_name>'",
            description="uncatalog a model service",
        )
    )

    statements.append(
        py.Forward(catalog + model + service + fr_om + quoted_string("path") + a_s + quoted_string("service_name"))(
            "catalog_add_model_service"
        )
    )
    help.append(
        help_dict_create(
            name="catalog Model servie",
            category="Model",
            command="catalog model service from '<path or github' as  '<service_name>'",
            description="catalog a model service from a path or github",
        )
    )

    statements.append(py.Forward(model + service + up + quoted_string("service_name"))("model_up"))
    help.append(
        help_dict_create(
            name="Model up",
            category="Model",
            command="model service up '<service_name>'",
            description="launch a model service",
        )
    )

    statements.append(py.Forward(model + service + down + quoted_string("service_name"))("model_down"))
    help.append(
        help_dict_create(
            name="Model down",
            category="Model",
            command="model service down '<service_name>'",
            description="bring down a model service",
        )
    )
