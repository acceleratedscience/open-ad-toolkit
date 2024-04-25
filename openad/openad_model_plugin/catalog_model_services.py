import pyparsing as py
import os
import json
import glob
from openad.helpers.output import output_text, output_table, output_warning, output_error, output_success
from openad.helpers.spinner import spinner
from openad.openad_model_plugin.services import ModelService, UserProvidedConfig
from typing import List, Dict, Tuple
from pandas  import DataFrame
from subprocess import run
import shlex
import shutil


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

    return output_table(DataFrame(ns), is_data=False, headers=["Cataloged Services"])


def model_service_status(cmd_pointer, parser):
    """This function catalogs a service"""
    # get list of directory names for the catalog models
    ns_status = []
    for service in Dispatcher.list():
        try:
            res = Dispatcher.get_short_status(service)
            status = ""
            if res.get("up"):
                service = f"<green>{service}</green>"
                status = "READY"
            elif res.get("url"):
                service = f"<yellow>{service}</yellow>"
                status = "PENDING"
            else:
                status = "DOWN"
            ns_status.append([service, status, res.get("url")])
        except Exception as e:
            # model service not cataloged or doesnt exist
            # output_warning(str(e) + ' (hint: catalog service must be installed and configured)')
            continue
    headers=["Service", "Status", "URL"]
    data = DataFrame(ns_status, columns=headers)
    if not ns_status:
        output_warning("No services available")
    else:
        output_table(ns_status, is_data=False, headers=headers)
    return data


def retrieve_model(from_path: str, to_path: str) -> Tuple[bool, str]:
    spinner.start("Retrieving model")
    # uses ssh or https
    if (from_path.startswith("git@") or from_path.startswith("https://")) and from_path.endswith(".git"):
        # test if git is available
        try:
            cmd = shlex.split("git --version")
            git_version = run(cmd, capture_output=True, text=True, check=True)
        except Exception as e:
            spinner.fail(f"git not installed or unreachable")
            spinner.stop()
            return False, "git not installed or unreachable"
        # attempt to download model using git ssh
        try:
            cmd = shlex.split(f"git clone {from_path} {to_path}")
            clone = run(cmd, capture_output=True, text=True)  # not running check=true
            assert clone.returncode == 0, clone.stderr
            spinner.succeed(f"successfully retrieved model {from_path}")
            spinner.stop()
            return True, ""
        except Exception as e:
            spinner.fail(f"error: {str(e)}")
            spinner.stop()
            return False, str(e)
    # uses local path
    elif os.path.exists(from_path):
        # attempt to copy model
        try:
            cmd = shlex.split(f"cp -r {from_path} {to_path}")
            cp = run(cmd, capture_output=True, text=True)
            assert cp.returncode == 0, cp.stderr
            spinner.succeed(f"successfully retrieved model {from_path}")
            spinner.stop()
            return True, ""
        except Exception as e:
            spinner.fail(f"failed to fetch path {from_path} >> {str(e)}")
            spinner.stop()
            return False, str(e)
    else:
        spinner.fail(f"invalid path {from_path}")
        spinner.stop()
        return False, f"invalid path {from_path}"

def catalog_add_model_service(cmd_pointer, parser):
    """Add model service repo to catalog"""
    service_name = parser.as_dict()["service_name"]
    path = parser.as_dict()["path"]
    # add service to api
    model_path = os.path.join(SERVICE_DEFINTION_PATH, service_name)
    
    is_model_path, _ = retrieve_model(path, model_path)
    
    if is_model_path is False: return
    with Dispatcher as service:
        # configure the sky yaml
        config = UserProvidedConfig(
            workdir=model_path,
            port=8090,
            setup="docker buildx build -f Dockerfile -t service .",
            run=f"docker run --rm --network host service",
            disk_size=100,
            )
        service.add_service(service_name, config)
    return output_success(f"service {service_name} added to catalog")


def service_up(cmd_pointer, parser) -> None:
    """This function synchronously starts a service"""
    service_name = parser.as_dict()["service_name"]
    spinner.start("Starting service")
    try:
        with Dispatcher as service:
            service.up(service_name)
        spinner.succeed(f"service ({service_name}) started")
    except Exception as e:
        spinner.fail(str(e))
    spinner.stop()
    # return output_success(f"service ({service_name}) started")


def start_service_shutdown(service_name):
    with Dispatcher as service:
        if service.status(service_name).get("url") or bool(service.status(service_name).get("up")):
            # shut down service
            service.down(service_name, force=True)
            # reinitialize service
            config = service.get_config(service_name)
            service.remove_service(service_name)
            service.add_service(service_name, config)
            output_warning(f"service {service_name} is terminating.. make take some time.")


def service_down(cmd_pointer, parser) -> None:
    """This function synchronously shuts down a service"""
    service_name = parser.as_dict()["service_name"]
    start_service_shutdown(service_name)


def uncatalog_model_service(cmd_pointer, parser):
    """This function removes a catalog from the ~/.openad_model_service directory"""
    service_name = parser.as_dict()["service_name"]
    if service_name not in Dispatcher.list():
        return output_error(f"service {service_name} not found in catalog")
    spinner.start(f"Removing service {service_name}")
    start_service_shutdown(service_name)
    with Dispatcher as service:
        try:
            shutil.rmtree(os.path.join(SERVICE_DEFINTION_PATH, service_name))
            service.remove_service(service_name)
            spinner.succeed(f"service {service_name} removed from catalog")
        except:
            spinner.fail(f"could not delete local service at {os.path.join(SERVICE_DEFINTION_PATH, service_name)}")
    spinner.stop()
    return


def get_service_endpoint(service_name) -> str | None:
    """gets the service endpoint for a given service, if endpoint is not available it returns None"""
    if service_name is None:
        # may in future return a default local service
        return None
    with Dispatcher as service:
        endpoint = json.loads(service.status(service_name)).get("url")
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
