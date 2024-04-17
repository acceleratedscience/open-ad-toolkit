import pyparsing as py
import os
import json
import glob
from openad.helpers.output import output_text, output_table, output_warning, output_error, output_success
from openad.openad_model_plugin.services import ModelService, UserProvidedConfig
from typing import List, Dict


SERVICE_DEFINTION_PATH = os.path.expanduser("~/.openad_model_services/")
SERVICES_PATH = os.path.expanduser("/definitions/services/")
if not os.path.exists(SERVICE_DEFINTION_PATH):
    os.makedirs(SERVICE_DEFINTION_PATH)

# initialize model services
model_service = ModelService()

model_service.download_model(name="gtsd4_gen", url="", model_dir=SERVICE_DEFINTION_PATH)
model_service.download_model(name="gtsd4_prop", url="", model_dir=SERVICE_DEFINTION_PATH)

### example of how to load services by namespace ###
# with model_service(TEST_PATH) as service:
#     print(service.list())

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


def model_service_status(cmd_pointer, parser) -> Dict:
    """Get a model service status"""
    service_name = parser.as_dict()["service_name"]
    with model_service(SERVICE_DEFINTION_PATH + service_name) as model:
        print(model.status(service_name))

def get_namespaces():
    list_of_namespaces = [
        os.path.basename(f.path) for f in os.scandir(SERVICE_DEFINTION_PATH) if f.is_dir()
    ]  # os.walk(SERVICE_DEFINTION_PATH)
    return list_of_namespaces
 
def list_cataloged_model_services(cmd_pointer, parser):
    """This function catalogs a service"""
    print("listing catalog services")
    namespaces = get_namespaces()
    ns_status = []
    print(f"{namespaces=}")
    for ns in namespaces:
        print("trying", ns)
        try:
            with model_service(SERVICE_DEFINTION_PATH + ns) as model:
                    ns_status.append({ns: model.get_short_status(ns)})
        except Exception as e:
            print('[Error]', e)
            continue  # service not cataloged or doesnt exist
    print(ns_status)
    # status = []
    # catalog = model_service.list()
    # for service in catalog:
    #     status.append({service: model_service.get_short_status(service)})
    # print(status)


def catalog_add_model_service(cmd_pointer, parser):
    """Add model service repo to catalog"""
    service_name = parser.as_dict()["service_name"]
    path = parser.as_dict()["path"]
    # add service to api
    model_path = SERVICE_DEFINTION_PATH + service_name
    print(f"{model_path=}")
    # model_service.save(model_path)  # add to catalog
    with model_service(model_path) as model:
        # configure the sky yaml
        config = UserProvidedConfig(
            workdir=model_path,
            run="python openad_model_property_service/property_service.py",
            # disk_size=100
            )
        model.add_service(service_name, config)
        print(f"service {service_name} added to catalog from {path}")


def service_up(cmd_pointer, parser) -> None:
    """This function synchronously starts a service"""
    service_name = parser.as_dict()["service_name"]
    print(SERVICE_DEFINTION_PATH + service_name)
    with model_service(SERVICE_DEFINTION_PATH + service_name) as model:
        model.up(service_name)
        print("service up")


def service_down(cmd_pointer, parser) -> None:
    """This function synchronously shuts down a service"""
    service_name = parser.as_dict()["service_name"]
    model_service.down(service_name)


def uncatalog_model_service(cmd_pointer, parser):
    """This function removes a catalog from the ~/.openad_model_service directory"""
    service_name = parser.as_dict()["service_name"]
    with model_service(SERVICE_DEFINTION_PATH + service_name) as model:
        try:
            model.down(service_name)
            model.remove_service(service_name)
        except Exception as e:
            print("dies")
            print(e)
        finally:
            print('cache: ', model.cache)
            print('services: ', model.list())
            model.remove_service(service_name)
            print("finally done")
            # model.save(SERVICE_DEFINTION_PATH + service_name)
    print(f"service {service_name} removed from catalog")


def get_service_endpoint(service_name) -> List:
    """gets the service endpoint for a given service, if endpoint is not available it returns None"""
    if service_name is None:
        "may in future return a default local service"
        return []
    endpoint = json.loads(model_service.status(service_name)).get("url")
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
