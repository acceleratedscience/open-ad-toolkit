import pyparsing as py
import os
import json
import glob
from openad.helpers.output import output_text, output_table, output_warning, output_error, output_success
from openad.helpers.spinner import spinner
from openad.openad_model_plugin.services import ModelService, UserProvidedConfig
from typing import List, Dict, Tuple
from pandas import DataFrame
from subprocess import run
import shlex
import shutil
from tabulate import tabulate
from tomlkit import parse
import time


SERVICE_DEFINTION_PATH = os.path.expanduser("~/.openad_model_services/")
SERVICES_PATH = "/definitions/services/"
if not os.path.exists(SERVICE_DEFINTION_PATH):
    os.makedirs(SERVICE_DEFINTION_PATH)


# this is the global object that should be used across openad and testing
Dispatcher = ModelService()
### example of how to use the dispatcher ###
# with Dispatcher() as service:
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
    service_files = glob.glob(reference + "/*.json", recursive=True)

    for file in service_files:
        with open(file, "r") as file_handle:
            try:
                jdoc = json.load(file_handle)
                service_list.append(jdoc)
            except Exception as e:
                output_error("invalid service json definition  " + file)
                output_error(e)
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
        else:
            services_path = SERVICE_DEFINTION_PATH + namespace + "/**" + SERVICES_PATH
            services_path = glob.glob(services_path, recursive=True)
            if len(services_path) > 0:
                services_path = services_path[0]
                service_list = get_service_defs(services_path)
        service_list_by_catalog[namespace] = service_list
    return service_list_by_catalog


def get_catalog_namespaces(cmd_pointer, parser) -> Dict:
    """Get a model service status"""
    ns = get_namespaces()

    return output_table(DataFrame(ns), headers=["Cataloged Services"], is_data=False)


def model_service_status(cmd_pointer, parser):
    """This function catalogs a service"""
    # get list of directory names for the catalog models
    models = {"Service": [], "Status": [], "Endpoint": []}
    with Dispatcher(update_status=True) as service:
        # get all the services then order by name and if url exists
        all_services: list = service.list()
        with_url: set = set(i for i in all_services if service.get_short_status(i).get("url"))
        without_url: set = set(all_services) - with_url
        order_services = sorted(list(with_url)) + sorted(list(without_url))
        # !important load services with update
        if all_services:  # proceed if any service available
            try:
                spinner.start("searching running services")
                # TODO: verify how much time or have a more robust method
                time.sleep(3)  # wait for service threads to ping endpoint
                for name in order_services:
                    res = service.get_short_status(name)
                    # set the status of the service
                    if res.get("up"):
                        status = "READY"
                    elif res.get("url"):
                        status = "PENDING"
                    else:
                        status = "DOWN"
                    models["Service"].append(name)
                    models["Status"].append(status)
                    if res.get("url"):
                        models["Endpoint"].append("http://" + res.get("url"))
                    else:
                        models["Endpoint"].append(res.get("url"))
            except Exception as e:
                # model service not cataloged or doesnt exist
                output_warning(str(e))
            finally:
                spinner.stop()
    return DataFrame(models)


def model_service_config(cmd_pointer, parser):
    """prints service resources"""
    service_name = parser.as_dict()["service_name"]
    with Dispatcher() as service:
        res = service.get_config_as_dict(service_name)
        config = {**res["template"]["service"], **res["template"]["resources"]}
        table_data = [[key, value] for key, value in config.items()]
        # print(tabulate(table_data, headers=["Resource", "value"], tablefmt="pretty"))
    return DataFrame(table_data, columns=["Resource", "value"])


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
            spinner.info(f"successfully retrieved model {from_path}")
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
            spinner.info(f"successfully retrieved model {from_path}")
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


def load_service_config(local_service_path: str) -> UserProvidedConfig:
    """loads service params from openad.cfg file"""
    cfg_map = {"port":int, "replicas":int, "cloud":str,"disk_size":int, "cpu":str, "memory":str, "accelerators":str, "setup":str, "run":str}
    if os.path.exists(os.path.join(local_service_path, "openad.cfg")):
        try:
            # open the document
            with open(os.path.join(local_service_path, "openad.cfg")) as f:
                parser = parse(f.read())
            conf = {}
            # check if [defaults] key exists if not ignore. allows for new fields in the future
            if "defaults" in parser.keys():
                parser = parser.get("defaults")
            # cast the values into a new dict
            for key, value in parser.items():
                key = key.lower()
                if value and key in cfg_map.keys():
                    conf[key] = cfg_map[key](value)  # cast the type to value
            # check if conf has any values
            if conf:
                # create a UserProvidedConfig with conf data
                spinner.info("found non defaults in openad.cfg")
                table_data = [[key, value] for key, value in conf.items()]
                print(tabulate(table_data, headers=["service spec", "value"], tablefmt="pretty"))
                return UserProvidedConfig(**conf, workdir=local_service_path)
            else:
                spinner.warn("error with (openad.cfg). Could not load user config. Loading defaults.")
        except Exception as e:
            output_error(str(e))
            spinner.warn("error with (openad.cfg). Could not load user config. Loading defaults.")
    # use default config
    return UserProvidedConfig(workdir=local_service_path)


def catalog_add_model_service(cmd_pointer, parser) -> bool:
    """Add model service repo to catalog"""
    service_name = parser.as_dict()["service_name"]
    remote_service = parser.as_dict()["path"]
    # check if service exists
    with Dispatcher() as service:
        if service_name in service.list():
            spinner.fail(f"service {service_name} already exists in catalog")
            output_error(f"service {service_name} already exists in catalog", return_val=False)
            return False
    # download model
    local_service_path = os.path.join(SERVICE_DEFINTION_PATH, service_name)
    is_local_service_path, _ = retrieve_model(remote_service, local_service_path)
    if is_local_service_path is False:
        spinner.fail(f"service {service_name} was unable to be added to check url or path")
        output_error(f"service {service_name} was unable to be added to check url or path", return_val=False)
        spinner.stop()
        return False
    # get any available configs from service
    config = load_service_config(local_service_path)
    # add the service
    with Dispatcher() as service:
        service.add_service(service_name, config)
        # spinner.succeed(f"service {service_name} added to catalog")
        output_success(f"service {service_name} added to catalog", return_val=False)
    return True


def uncatalog_model_service(cmd_pointer, parser):
    """This function removes a catalog from the ~/.openad_model_service directory"""
    service_name = parser.as_dict()["service_name"]
    with Dispatcher() as service:
        # check if service exists
        if service_name not in service.list():
            return output_error(f"service {service_name} not found in catalog", return_val=False)
            return False
        # stop running service
        start_service_shutdown(service_name)
        # remove local files for service
        if os.path.exists(os.path.join(SERVICE_DEFINTION_PATH, service_name)):
            shutil.rmtree(os.path.join(SERVICE_DEFINTION_PATH, service_name))
        # remove service from cache
    with Dispatcher() as service:  # initialize fresh load
        try:
            service.remove_service(service_name)
            spinner.succeed(f"service {service_name} removed from catalog")
        except Exception as e:
            if "No such file or directory" in str(e):
                spinner.warn("service doesnt exist but trying to remove from list. config file was already deleted")
                # TODO: make more robust error handling
                path = os.path.join(os.path.expanduser("~/.servicing"), f"{service_name}_service.yaml")
                open(path).close()  # create file
                service.remove_service(service_name)
            else:
                spinner.fail(f"failed to remove service: {str(e)}")
                # output_error(f"failed to remove service: {str(e)}", return_val=False)
                return False
        output_success(f"service {service_name} removed from catalog", return_val=False)
    return True


def service_up(cmd_pointer, parser) -> None:
    """This function synchronously starts a service"""
    service_name = parser.as_dict()["service_name"]
    # spinner.start("Starting service")
    try:
        with Dispatcher() as service:
            service.up(service_name, skip_prompt=True)
        # spinner.succeed(f"service ({service_name}) started")
    except Exception as e:
        output_error("Service was unable to be started:\n" + str(e), return_val=False)
        return False
    output_success(f"Service {service_name} is Starting.. may take some time.", return_val=False)
    return True
    # spinner.stop()
    # return output_success(f"service ({service_name}) started")


def service_up_endpoint(cmd_pointer, parser) -> None:
    endpoint = parser.as_dict()["endpoint"]
    output_error("Not yet implemented")


def local_service_up(cmd_pointer, parser) -> None:
    service_name = parser.as_dict()["service_name"]
    output_error("Not yet implemented")


def start_service_shutdown(service_name):
    with Dispatcher() as service:
        if service.status(service_name).get("url") or bool(service.status(service_name).get("up")):
            # shut down service
            service.down(service_name, skip_prompt=True)
            # reinitialize service
            config = service.get_user_provided_config(service_name)
            service.remove_service(service_name)
            service.add_service(service_name, config)
            spinner.warn(f"service {service_name} is terminating.. may take some time.")
            return True
        else:
            # output_error(
            #    f"service {service_name} was not able to terminate, please check error sky pilot to determine status and force shutdown",
            #    return_val=False,
            # )
            return False


def service_down(cmd_pointer, parser) -> None:
    """This function synchronously shuts down a service"""
    is_success = False
    try:
        service_name = parser.as_dict()["service_name"]
        spinner.start(f"terminating {service_name} service")
        if not start_service_shutdown(service_name):
            spinner.info(f"service {service_name} is not up")
            # output_warning(f"service {service_name} is not up")
            is_success = True
    except Exception as e:
        output_error(str(e))
    finally:
        spinner.stop()
    return is_success


def get_service_endpoint(service_name) -> str | None:
    """gets the service endpoint for a given service, if endpoint is not available it returns None"""
    if service_name is None:
        # may in future return a default local service
        return None
    with Dispatcher() as service:

        # endpoint = json.loads(service.status(service_name)).get("url")
        endpoint = service.status(service_name)["url"]

    return endpoint


def service_catalog_grammar(statements: list, help: list):
    """This function creates the required grammar for managing cataloging services and model up or down"""
    catalog = py.CaselessKeyword("catalog")
    uncatalog = py.CaselessKeyword("uncatalog")
    model = py.CaselessKeyword("model")
    up = py.CaselessKeyword("up")
    local = py.CaselessKeyword("local")
    down = py.CaselessKeyword("down")
    service = py.CaselessKeyword("service")
    status = py.CaselessKeyword("status")
    fr_om = py.CaselessKeyword("from")
    _list = py.CaselessKeyword("list")
    path = py.CaselessKeyword("path")
    quoted_string = py.QuotedString("'", escQuote="\\")
    a_s = py.CaselessKeyword("as")
    config = py.CaselessKeyword("config")
    remote = py.CaselessKeyword("remote")

    statements.append(py.Forward(model + service + status)("model_service_status"))
    help.append(
        help_dict_create(
            name="model service status",
            category="Model",
            command="model service status",
            description="get the status of currently cataloged services",
        )
    )

    statements.append(py.Forward(model + service + up + remote + quoted_string("endpoint"))("service_up_endpoint"))
    help.append(
        help_dict_create(
            name="model service up remote",
            category="Model",
            command="model service up remote <endpoint>",
            description="connect to a remote model endpoint",
        )
    )

    statements.append(py.Forward(model + service + config + quoted_string("service_name"))("model_service_config"))
    help.append(
        help_dict_create(
            name="model service config",
            category="Model",
            command="model service config '<service_name>'",
            description="get the config of a service",
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
            name="catalog Model service",
            category="Model",
            command="catalog model service from '<path or github>' as  '<service_name>'",
            description="catalog a model service from a path or github",
        )
    )

    statements.append(py.Forward(model + service + up + quoted_string("service_name"))("service_up"))
    help.append(
        help_dict_create(
            name="Model up",
            category="Model",
            command="model service up '<service_name>'",
            description="launch a model service",
        )
    )

    statements.append(py.Forward(model + service + local + up + quoted_string("service_name"))("local_service_up"))
    help.append(
        help_dict_create(
            name="Model local up",
            category="Model",
            command="model service local up '<service_name>'",
            description="launch a model service locally",
        )
    )

    statements.append(py.Forward(model + service + down + quoted_string("service_name"))("service_down"))
    help.append(
        help_dict_create(
            name="Model down",
            category="Model",
            command="model service down '<service_name>'",
            description="bring down a model service",
        )
    )
