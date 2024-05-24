import pyparsing as py
import os
import sys
import json
import glob
from openad.helpers.output import (
    output_text,
    output_table,
    output_warning,
    output_error,
    output_success,
)
from openad.helpers.spinner import spinner
from openad.app.global_var_lib import GLOBAL_SETTINGS
from openad.openad_model_plugin.services import ModelService, UserProvidedConfig, ServiceFetchError
from typing import List, Dict, Tuple
from pandas import DataFrame
from subprocess import run
import shlex
import shutil
from tabulate import tabulate
from tomlkit import parse
import time
from openad.openad_model_plugin.utils import get_logger, bcolors
from functools import cache, lru_cache


logger = get_logger(__name__, color=bcolors.OKCYAN + bcolors.UNDERLINE)


DISPATCHER_SERVICE_PATH = os.path.expanduser("~/.servicing/")
SERVICE_MODEL_PATH = os.path.expanduser("~/.openad_model_services/")
SERVICES_PATH = "/definitions/services/"
if not os.path.exists(SERVICE_MODEL_PATH):
    os.makedirs(SERVICE_MODEL_PATH)


# this is the global object that should be used across openad and testing
logger.debug("initializing global Model Service.")
Dispatcher = ModelService(location=DISPATCHER_SERVICE_PATH, update_status=True, skip_sky_validation=True)
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
        os.path.basename(f.path) for f in os.scandir(SERVICE_MODEL_PATH) if f.is_dir()
    ]  # os.walk(SERVICE_MODEL_PATH)
    logger.debug(f"finding namespaces | {list_of_namespaces=}")
    return list_of_namespaces


@lru_cache(maxsize=16)
def get_local_service_defs(reference: str) -> list:
    """pulls the list of available service definitions. caches first result"""
    logger.debug(f"searching defs in {reference}")
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


def load_service_cache() -> Dict[str, dict] | None:
    """load latest cache"""
    # TODO: implement load at beginning. need to make custom cache
    try:
        with open(os.path.join(DISPATCHER_SERVICE_PATH, "cached_service_defs.json"), "w+") as f_cached:
            logger.debug("loading service defs cache")
            return json.load(f_cached)
    except Exception as e:
        # catch all. if it fails should not stop us from proceeding
        logger.error(f"could not load service defs cache: {str(e)}")
        return None


def save_service_cache(service_definitions):
    """save latest cache"""
    try:
        with open(os.path.join(DISPATCHER_SERVICE_PATH, "cached_service_defs.json"), "w+") as f_cached:
            logger.debug("saving cache to file")
            json.dump(service_definitions, f_cached)
    except Exception as e:
        # catch all. if it fails should not stop us from proceeding
        logger.error(f"could not save service defs cache to file: {str(e)}")
        pass


def get_cataloged_service_defs() -> Dict[str, dict]:
    """Returns a dictionary of cataloged services definitions"""
    logger.debug("checking available service definitions")
    service_definitions = dict()
    # get local namespace service definitions
    list_of_namespaces = [os.path.basename(f.path) for f in os.scandir(SERVICE_MODEL_PATH) if f.is_dir()]
    # iterate over local service definitions
    for namespace in list_of_namespaces:
        service_list = []
        services_path = SERVICE_MODEL_PATH + namespace + SERVICES_PATH
        if os.path.exists(services_path):
            service_list = get_local_service_defs(services_path)
        else:
            services_path = SERVICE_MODEL_PATH + namespace + "/**" + SERVICES_PATH
            services_path = glob.glob(services_path, recursive=True)
            if len(services_path) > 0:
                services_path = services_path[0]
                service_list = get_local_service_defs(services_path)
        if service_list:
            logger.debug(f"adding local defs for | {namespace=}")
            service_definitions[namespace] = service_list
    # iterate over remote service definitions
    with Dispatcher() as service:
        dispatcher_services = service.list()
        # iterate over keys not used before
        for name in set(dispatcher_services) - set(list_of_namespaces):
            remote_definitions = service.get_remote_service_definitions(name)
            if remote_definitions:
                logger.debug(f"adding remote service defs for | {name=}")
                service_definitions[name] = remote_definitions
    return service_definitions


def get_catalog_namespaces(cmd_pointer, parser) -> Dict:
    """Get a local model catalog"""
    ns = get_namespaces()
    return output_table(DataFrame(ns), headers=["Cataloged Services"], is_data=False)


def model_service_status(cmd_pointer, parser):
    """get all services status"""
    logger.debug("listing model status")
    # get list of directory names for the catalog models
    models = {"Service": [], "Status": [], "Endpoint": [], "Type": []}
    with Dispatcher(update_status=True) as service:
        # get all the services then order by name and if url exists
        all_services: list = service.list()
        # !important load services with update
        if all_services:  # proceed if any service available
            try:
                spinner.start("searching running services")
                # TODO: verify how much time or have a more robust method
                time.sleep(2)  # wait for service threads to ping endpoint
                for name in all_services:
                    res = service.get_short_status(name)
                    # set the status of the service
                    if res.get("up"):
                        status = "READY"
                    elif res.get("url") and not res.get("is_remote"):
                        status = "PENDING"
                    elif res.get("is_remote") and res.get("url"):
                        status = "UNREACHABLE"
                    else:
                        status = "DOWN"
                    if res.get("is_remote"):
                        models["Type"].append("remote")
                    else:
                        models["Type"].append("local")
                    models["Service"].append(name)
                    models["Status"].append(status)
                    models["Endpoint"].append(res.get("url"))
            except Exception as e:
                # model service not cataloged or doesnt exist
                output_warning(f"Error getting status: {str(e)}")
            finally:
                spinner.stop()
    df = DataFrame(models)
    return df.sort_values(by=["Status", "Service"], ascending=[False, True])


def model_service_config(cmd_pointer, parser):
    """prints service resources"""
    logger.debug("listing service config")
    service_name = parser.as_dict()["service_name"]
    with Dispatcher() as service:
        res = service.get_config_as_dict(service_name)
        config = {**res["template"]["service"], **res["template"]["resources"]}
        table_data = [[key, value] for key, value in config.items()]
        # print(tabulate(table_data, headers=["Resource", "value"], tablefmt="pretty"))
    return DataFrame(table_data, columns=["Resource", "value"])


def retrieve_model(from_path: str, to_path: str) -> Tuple[bool, str]:
    logger.debug("retrieving service model")
    spinner.start("Retrieving model")
    # uses ssh or https
    if (from_path.startswith("git@") or from_path.startswith("https://")) and from_path.endswith(".git"):
        # test if git is available
        try:
            cmd = shlex.split("git --version")
            run(cmd, capture_output=True, text=True, check=True)
        except Exception:
            spinner.fail("git not installed or unreachable")
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
    logger.debug(f"get local service configuration | {local_service_path=}")
    cfg_map = {
        "port": int,
        "replicas": int,
        "cloud": str,
        "disk_size": int,
        "cpu": str,
        "memory": str,
        "accelerators": str,
        "setup": str,
        "run": str,
    }
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
                print(tabulate(table_data, headers=["Resource", "value"], tablefmt="pretty"))
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
    service_path = os.path.expanduser(parser.as_dict()["path"])
    logger.debug(f"catalog model service | {service_name=} {service_path=}")
    if "remote" in parser:
        return service_up_endpoint(cmd_pointer, parser)
    # check if service exists
    with Dispatcher() as service:
        if service_name in service.list():
            spinner.fail(f"service {service_name} already exists in catalog")
            output_error(f"service {service_name} already exists in catalog", return_val=False)
            return False
    # download model
    local_service_path = os.path.join(SERVICE_MODEL_PATH, service_name)
    is_local_service_path, _ = retrieve_model(service_path, local_service_path)
    if is_local_service_path is False:
        spinner.fail(f"service {service_name} was unable to be added to check url or path")
        output_error(
            f"service {service_name} was unable to be added to check url or path",
            return_val=False,
        )
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


def uncatalog_model_service(cmd_pointer, parser) -> bool:
    """This function removes a catalog from the ~/.openad_model_service directory"""
    service_name = parser.as_dict()["service_name"]
    logger.debug(f"uncatalog model service | {service_name=}")
    with Dispatcher() as service:
        # check if service exists
        if service_name not in service.list():
            output_error(f"service {service_name} not found in catalog", return_val=False)
            return False
        # stop running service
        start_service_shutdown(service_name)
        # remove local files for service
        if os.path.exists(os.path.join(SERVICE_MODEL_PATH, service_name)):
            shutil.rmtree(os.path.join(SERVICE_MODEL_PATH, service_name))
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


def service_up(cmd_pointer, parser) -> bool:
    """This function synchronously starts a service"""
    gpu_disable = "no_gpu" in parser.as_dict()  # boolean flag to disable gpu
    service_name = parser.as_dict()["service_name"]
    logger.debug(f"start service | {service_name=} {gpu_disable=}")
    # spinner.start("Starting service")
    output_success("Deploying Service. Please Wait.....", return_val=False)
    try:
        with Dispatcher() as service:
            service.up(service_name, skip_prompt=True, gpu_disable=gpu_disable)

            # spinner.succeed(f"service ({service_name}) started")
            output_success(f"Service {service_name} is Starting.. may take some time.")
            return True
    except Exception as e:
        output_error("Service was unable to be started:\n" + str(e), return_val=False)
        return False


def service_up_endpoint(cmd_pointer, parser) -> bool:
    service_name = parser.as_dict()["service_name"]
    endpoint = parser.as_dict()["path"]
    logger.debug(f"add as remote service | {service_name=} {endpoint=}")
    with Dispatcher() as service:
        if service_name in service.list():
            spinner.fail(f"service {service_name} already exists")
            return False
        # load remote endpoint to config custom field
        config = json.dumps(
            {
                "remote_service": True,
                "remote_endpoint": endpoint,
                "remote_status": False,
            }
        )
        service.add_service(service_name, UserProvidedConfig(data=config))
    spinner.succeed(f"Remote service '{service_name}' added!")
    return True


def local_service_up(cmd_pointer, parser) -> None:
    service_name = parser.as_dict()["service_name"]
    logger.debug(f"start service locally | {service_name=}")
    output_error(f" {service_name} Not yet implemented")


def start_service_shutdown(service_name):
    logger.debug(f"prepare service shutdown | {service_name=}")
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
        logger.debug(f"attempt to stop service | {service_name=}")
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
        endpoint = service.get_url(service_name)
    logger.debug(f"get service endpoint | {service_name=} {endpoint=}")
    return endpoint


def service_catalog_grammar(statements: list, help: list):
    """This function creates the required grammar for managing cataloging services and model up or down"""
    logger.debug("catalog model service grammer")
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

    statements.append(
        py.Forward(model + service + config + (quoted_string | py.Word(py.alphanums + "_"))("service_name"))(
            "model_service_config"
        )
    )
    help.append(
        help_dict_create(
            name="model service config",
            category="Model",
            command="model service config '<service_name>'|<service_name>",
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
        py.Forward(uncatalog + model + service + (quoted_string | py.Word(py.alphanums + "_"))("service_name"))(
            "uncatalog_model_service"
        )
    )
    help.append(
        help_dict_create(
            name="uncatalog model service",
            category="Model",
            command="uncatalog model service '<service_name>'|<service_name>",
            description="uncatalog a model service \n\n Example: \n<cmd>uncatalog model service 'gen'</cmd>",
        )
    )

    statements.append(
        py.Forward(
            catalog
            + model
            + service
            + fr_om
            + py.Optional(remote("remote"))
            + quoted_string("path")
            + a_s
            + (quoted_string | py.Word(py.alphanums + "_"))("service_name")
        )("catalog_add_model_service")
    )
    help.append(
        help_dict_create(
            name="catalog Model service",
            category="Model",
            command="catalog model service from (remote) '<path or github>' as  '<service_name>'|<service_name>",
            description="""catalog a model service from a path or github or remotely from an existing OpenAD service.

Example:

-<cmd>catalog model service from 'git@github.com:acceleratedscience/generation_inference_service.git' as 'gen'</cmd>

or to catalog a remote service shared with you:
-<cmd>catalog model service from remote 'http://54.235.3.243:30001' as gen</cmd>""",
        )
    )

    statements.append(
        py.Forward(
            model
            + service
            + up
            + (quoted_string | py.Word(py.alphanums + "_"))("service_name")
            + py.Optional(py.CaselessKeyword("NO_GPU")("no_gpu"))
        )("service_up")
    )
    help.append(
        help_dict_create(
            name="Model up",
            category="Model",
            command="model service up '<service_name>'|<service_name> [no_gpu]}",
            description="""launches a cataloged model service.
If you do not want to launch a service with GPU you should specify <cmd>no_gpu</cmd> at the end of the command.
Examples:

-<cmd>model service up gen</cmd>

-<cmd>model service up 'gen'</cmd>

-<cmd>model service up gen no_gpu</cmd>""",
        )
    )

    statements.append(
        py.Forward(
            model
            + service
            + local
            + up
            + (quoted_string | py.Word(py.alphanums + "_"))("service_name")
            + py.Optional(py.CaselessKeyword("NO_GPU")("no_gpu"))
        )("local_service_up")
    )
    help.append(
        help_dict_create(
            name="Model local up",
            category="Model",
            command="model service local up '<service_name>'|<service_name> ",
            description="""launch a model service locally.

            Example:
              <cmd> model service local up gen</cmd>

             """,
        )
    )

    statements.append(
        py.Forward(model + service + down + (quoted_string | py.Word(py.alphanums + "_"))("service_name"))(
            "service_down"
        )
    )
    help.append(
        help_dict_create(
            name="Model down",
            category="Model",
            command="model service down '<service_name>'|<service_name>",
            description="Bring down a model service  \n Examples: \n\n<cmd>model service down gen</cmd> \n\n<cmd>model service down 'gen'</cmd> ",
        )
    )
