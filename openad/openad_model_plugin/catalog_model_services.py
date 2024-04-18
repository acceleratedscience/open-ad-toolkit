import pyparsing as py
import os
import json
import glob
from openad.helpers.output import output_text, output_table, output_warning, output_error, output_success
from openad.openad_model_plugin.services import ModelService, UserProvidedConfig
from typing import List, Dict
from pandas  import DataFrame


SERVICE_DEFINTION_PATH = os.path.expanduser("~/.openad_model_services/")
SERVICES_PATH = os.path.expanduser("/definitions/services/")
if not os.path.exists(SERVICE_DEFINTION_PATH):
    os.makedirs(SERVICE_DEFINTION_PATH)


Dispatcher = ModelService()
### example of how to use the dispatcher ###
# with Dispatcher as service:
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


def get_namespaces():
    list_of_namespaces = [
        os.path.basename(f.path) for f in os.scandir(SERVICE_DEFINTION_PATH) if f.is_dir()
    ]  # os.walk(SERVICE_DEFINTION_PATH)
    return list_of_namespaces


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

    return output_table(DataFrame(ns), is_data=False, headers=["Services"])


def model_service_status(cmd_pointer, parser):
    """This function catalogs a service"""
    # get list of directory names for the catalog models
    namespaces = get_namespaces()
    ns_status = []
    for ns in namespaces:
        try:
            res = Dispatcher.get_short_status(ns)
            status = ""
            if res.get("up"):
                ns = f"<green>{ns}</green>"
                status = "READY"
            elif res.get("url"):
                ns = f"<yellow>{ns}</yellow>"
                status = "PENDING"
            else:
                status = "DOWN"
            ns_status.append([ns, status, res.get("url")])
        except Exception as e:
            # model service not cataloged or doesnt exist
            # output_warning(str(e) + ' (hint: catalog service must be installed and configured)')
            continue
    return output_table(ns_status, is_data=False, headers=["Service", "Status", "URL"])


def catalog_add_model_service(cmd_pointer, parser):
    """Add model service repo to catalog"""
    service_name = parser.as_dict()["service_name"]
    path = parser.as_dict()["path"]
    # add service to api
    model_path = SERVICE_DEFINTION_PATH + service_name
    with Dispatcher as service:
        # configure the sky yaml
        config = UserProvidedConfig(
            workdir=model_path,
            port=8080,
            setup="docker buildx build -f Dockerfile -t service .",
            run="docker run --rm \
                    -p 8080:8080 \
                    service",
            disk_size=100,
            )
        service.add_service(service_name, config)
    return output_success(f"service {service_name} added to catalog")


def service_up(cmd_pointer, parser) -> None:
    """This function synchronously starts a service"""
    service_name = parser.as_dict()["service_name"]
    with Dispatcher as service:
        service.up(service_name)
    return output_success(f"service ({service_name}) started")


def service_down(cmd_pointer, parser) -> None:
    """This function synchronously shuts down a service"""
    service_name = parser.as_dict()["service_name"]
    with Dispatcher as service:
        service.down(service_name)
    return output_success(f"service {service_name} terminating..")


def uncatalog_model_service(cmd_pointer, parser):
    """This function removes a catalog from the ~/.openad_model_service directory"""
    service_name = parser.as_dict()["service_name"]
    with Dispatcher as service:
        try:
            service.down(service_name)
            output_warning(f"service {service_name} terminating..")
        except Exception as e:
            pass
        finally:
            service.remove_service(service_name)
    return output_success(f"service {service_name} removed from catalog")


def get_service_endpoint(service_name) -> List:
    """gets the service endpoint for a given service, if endpoint is not available it returns None"""
    if service_name is None:
        "may in future return a default local service"
        return []
    endpoint = json.loads(Dispatcher.status(service_name)).get("url")
    return output_success(f"{endpoint}")


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
