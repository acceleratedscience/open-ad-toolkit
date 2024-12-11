import glob
import json
import os
import shlex
import shutil
import sys
import time
from functools import lru_cache
from subprocess import run
from typing import Dict, Tuple
import pandas as pd

import pyparsing as py
from openad.app.global_var_lib import GLOBAL_SETTINGS
from openad.core.help import help_dict_create
from openad.helpers.output import output_error, output_success, output_table, output_text, output_warning
from openad.helpers.spinner import spinner
from openad.openad_model_plugin.auth_services import (
    load_lookup_table,
    remove_auth_group,
    remove_service_group,
    update_lookup_table,
    get_service_api_key,
    hide_api_keys,
)
from openad.openad_model_plugin.config import DISPATCHER_SERVICE_PATH, SERVICE_MODEL_PATH, SERVICES_PATH
from openad.openad_model_plugin.services import ModelService, UserProvidedConfig
from openad.openad_model_plugin.utils import bcolors, get_logger
from pandas import DataFrame
from tabulate import tabulate
from tomlkit import parse

logger = get_logger(__name__, color=bcolors.OKCYAN + bcolors.UNDERLINE)


# this is the global object that should be used across openad and testing
logger.debug("initializing global Model Service.")
Dispatcher = ModelService(location=DISPATCHER_SERVICE_PATH, update_status=True, skip_sky_validation=True)
### example of how to use the dispatcher ###
# with Dispatcher() as service:
#     print(service.list())


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
    # list_of_namespaces = [os.path.basename(f.path) for f in os.scandir(SERVICE_MODEL_PATH) if f.is_dir()]
    # list_of_namespaces = []
    # # iterate over local service definitions
    # for namespace in list_of_namespaces:
    #     service_list = []
    #     services_path = SERVICE_MODEL_PATH + namespace + SERVICES_PATH
    #     if os.path.exists(services_path):
    #         service_list = get_local_service_defs(services_path)
    #     else:
    #         services_path = SERVICE_MODEL_PATH + namespace + "/**" + SERVICES_PATH
    #         services_path = glob.glob(services_path, recursive=True)
    #         if len(services_path) > 0:
    #             services_path = services_path[0]
    #             service_list = get_local_service_defs(services_path)
    #     if service_list:
    #         logger.debug(f"adding local defs for | {namespace=}")
    #         service_definitions[namespace] = service_list
    # iterate over remote service definitions
    with Dispatcher() as service:
        dispatcher_services = service.list()
        # iterate over keys not used before
        for name in set(dispatcher_services):  # - set(list_of_namespaces):
            remote_definitions = service.get_remote_service_definitions(name)
            if remote_definitions:
                logger.debug(f"adding remote service defs for | {name=}")
                service_definitions[name] = remote_definitions
            else:
                logger.warning(f"remote service defs not found, sevice not available | {name=}")
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
    models = {"Service": [], "Status": [], "Endpoint": [], "Host": [], "API expires": []}
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
                    if res.get("message"):
                        # an overwite if something occured
                        status = res.get("message")
                    elif res.get("up"):
                        status = "Ready"
                    elif res.get("url") and not res.get("is_remote"):
                        status = "Pending"
                    elif res.get("is_remote") and res.get("url"):
                        status = "Unreachable"
                    else:
                        status = "DOWN"
                    if res.get("is_remote"):
                        models["Host"].append("remote")
                        proxy_info: dict = res.get("jwt_info")
                        models["API expires"].append(proxy_info.get("exp_formatted", "No Info"))
                    else:
                        models["Host"].append("local")
                        models["API expires"].append("")
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
    # load service status details
    with Dispatcher() as service:
        res = service.get_config_as_dict(service_name)
        config = {**res["template"]["service"], **res["template"]["resources"]}
        table_data = [[key, value] for key, value in config.items()]
    # add authentication group details
    auth_lookup_table = load_lookup_table()
    table_data.insert(0, ["authentication group", auth_lookup_table["service_table"].get(service_name, "None")])
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
                spinner.stop()
                table_data = [[key, value] for key, value in conf.items()]
                print(tabulate(table_data, headers=["Resource", "value"], tablefmt="pretty"))
                return UserProvidedConfig(**conf, workdir=local_service_path, data=json.dumps({}))
            else:
                spinner.warn("error with (openad.cfg). Could not load user config. Loading defaults.")
        except Exception as e:
            output_error(str(e))
            spinner.warn("error with (openad.cfg). Could not load user config. Loading defaults.")
            spinner.stop()
    # use default config
    return UserProvidedConfig(workdir=local_service_path, data=json.dumps({}))


def add_remote_service_from_endpoint(cmd_pointer, parser) -> bool:
    service_name = parser.as_dict()["service_name"]
    endpoint = parser.as_dict()["path"]
    logger.debug(f"add as remote service | {service_name=} {endpoint=}")
    with Dispatcher() as service:
        if service_name in service.list():
            spinner.fail(f"service {service_name} already exists")
            return False
        # load remote endpoint to config custom field
        if "params" in parser:
            params = {k: v for k, v in parser.as_dict().get("params")}
            logger.debug(f"user added params: {params}")
        else:
            params = {}
        config = json.dumps(
            {
                "remote_service": True,
                "remote_endpoint": endpoint,
                "remote_status": False,
                "params": params,  # header values for request
            }
        )
        service.add_service(service_name, UserProvidedConfig(data=config))
    spinner.succeed(f"Remote service '{service_name}' added!")
    return True


def catalog_add_model_service(cmd_pointer, parser) -> bool:
    """Add model service repo to catalog"""

    service_name = parser.as_dict()["service_name"]
    service_path = os.path.expanduser(parser.as_dict()["path"])
    logger.debug(f"catalog model service | {service_name=} {service_path=}")
    params = {}
    if "params" in parser.as_dict():
        for i in parser.as_dict()["params"]:
            params[i[0]] = i[1]

    if "auth_group" in params.keys() and "authorization" in params.keys():
        output_error(
            f"It is not permitted to define auth_group and authroization in the same catalog statement ",
            return_val=False,
        )
        return False

    if "auth_group" in params.keys():
        auth_group = params["auth_group"]
        lookup_table = load_lookup_table()
        if auth_group not in lookup_table["auth_table"]:
            output_error(
                f"auth_group {auth_group} not in auth table, please add the authgroup and recatalog the service {service_name} ",
                return_val=False,
            )
            return False
    else:
        auth_group = None

    if "remote" in parser:
        # run this code and exit
        output = add_remote_service_from_endpoint(cmd_pointer, parser)
        if auth_group is not None:
            updated_lookup_table = update_lookup_table(auth_group=auth_group, service=service_name)
        output_success(f"Service {service_name} added to catalog for remote service {service_path}", return_val=False)
        return output
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
        output_success(f"Service {service_name} added to catalog", return_val=False)
    # If auth group in parameters apply authgroup

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
        output_success(f"Service {service_name} removed from catalog", return_val=False)
    # remove service from authentication lookup table
    if get_service_api_key(service_name):
        remove_service_group(service_name)

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
            spinner.stop()
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


def get_service_requester(service_name) -> str | None:
    """gets the service request params for a given service, if endpoint is not available it returns None"""
    if service_name is None:
        # may in future return a default local service
        return None
    with Dispatcher() as service:
        status = service.get_short_status(service_name)
        endpoint = service.get_url(service_name)
        return {"func": service.service_request, "status": status, "endpoint": endpoint}


def add_service_auth_group(cmd_pointer, parser):
    """Create an authentication group"""
    auth_group = parser.as_dict()["auth_group"]
    api_key = parser.as_dict()["api_key"]
    logger.debug(f"adding auth group | {auth_group=} {api_key=}")
    lookup_table = load_lookup_table()
    if auth_group in lookup_table["auth_table"]:
        return output_error(f"authentication group '{auth_group}' already exists")
    updated_lookup_table = update_lookup_table(auth_group=auth_group, api_key=api_key)
    output_success(f"successfully added authentication group '{auth_group}'")
    hide_api_keys(updated_lookup_table)
    return DataFrame(updated_lookup_table["auth_table"].items(), columns=["auth group", "api key"])


def remove_service_auth_group(cmd_pointer, parser):
    """remove an authentication group"""
    auth_group = parser.as_dict()["auth_group"]
    logger.debug(f"removing auth group | {auth_group=}")
    lookup_table = load_lookup_table()
    if auth_group not in lookup_table["auth_table"]:
        return output_error(f"authentication group '{auth_group}' does not exists")
    updated_lookup_table = remove_auth_group(auth_group)
    output_success(f"removed authentication group '{auth_group}'")
    hide_api_keys(updated_lookup_table)
    return DataFrame(updated_lookup_table["auth_table"].items(), columns=["auth group", "api key"])


def attach_service_auth_group(cmd_pointer, parser):
    """add a model service to an authentication group"""
    service_name = parser.as_dict()["service_name"]
    auth_group = parser.as_dict()["auth_group"]
    logger.debug(f"attaching auth group to service | {service_name=} {auth_group=}")
    lookup_table = load_lookup_table()
    # connect mapping to service from auth group
    with Dispatcher() as dispatch:
        models = dispatch.list()
        if service_name not in models:
            return output_error(f"service '{service_name}' does not exist")
        if auth_group not in lookup_table["auth_table"]:
            return output_error(f"auth group '{auth_group}' does not exist")
    # add auth to service
    updated_lookup_table = update_lookup_table(auth_group=auth_group, service=service_name)
    hide_api_keys(updated_lookup_table)
    return DataFrame(updated_lookup_table["service_table"].items(), columns=["service", "auth group"])


def detach_service_auth_group(cmd_pointer, parser):
    """remove a model service from an authentication group"""
    service_name = parser.as_dict()["service_name"]
    logger.debug(f"detaching auth group from service | {service_name=}")
    lookup_table = load_lookup_table()
    if service_name not in lookup_table["service_table"]:
        return output_error(f"service '{service_name}' does not have an authentication group")
    updated_lookup_table = remove_service_group(service_name)
    hide_api_keys(updated_lookup_table)
    return DataFrame(updated_lookup_table["service_table"].items(), columns=["service", "auth group"])


def list_auth_services(cmd_pointer, parser):
    """list authentication groups and services that use it"""
    # Extracting the data from the dictionary
    lookup_table = load_lookup_table(hide_api=True)
    services = []
    auth_groups = []
    apis = []
    # Extract services and their corresponding auth groups
    for service, auth_group in lookup_table["service_table"].items():
        services.append(service)
        auth_groups.append(auth_group)
        apis.append(lookup_table["auth_table"].get(auth_group))
    # Add auth groups from auth_table that are not in service_table
    for auth_group, api in lookup_table["auth_table"].items():
        if auth_group not in auth_groups:
            services.append(None)
            auth_groups.append(auth_group)
            apis.append(api)
    # Creating the DataFrame
    return DataFrame({"service": services, "auth group": auth_groups, "api key": apis})


def get_model_service_result(cmd_pointer, parser):
    # with Dispatcher as servicer:
    #    service_status = servicer.get_short_status(parser.to_dict()["service_name"].lower())
    try:
        # response = Dispatcher.service_request(
        #     name=service_name, method="POST", timeout=None, verify=not service_status.get("is_remote"), _json=a_request
        # )
        a_request = {"url": parser.as_dict()["request_id"], "service_type": "get_result"}
        response = Dispatcher.service_request(
            name=parser.as_dict()["service_name"].lower(), method="POST", timeout=None, verify=False, _json=a_request
        )
        # response = requests.post(Endpoint + "/service", json=a_request, headers=headers, verify=False)
    except Exception as e:
        output_error(str(e))
        return output_error("Error: \n Server not reachable at ")

    try:
        response_result = response.json()
        try:
            if isinstance(response_result, str):
                response_result = json.loads(response_result)
            if isinstance(response_result, dict):
                if "warning" in response_result:
                    return output_warning(response_result["warning"]["reason"])
                elif "error" in response_result:
                    run_error = "Request Error:\n"

                    for key, value in response_result["error"].items():
                        value = str(value).replace("<", "`<")
                        value = str(value).replace(">", ">`")
                        run_error = run_error + f"- <cmd>{key}</cmd> : {value}\n  "
                    return output_error(run_error)
                if "detail" in response_result:
                    return output_warning(response_result["detail"])

            result = pd.DataFrame(response_result)
            if "save_as" in parser:
                results_file = str(parser["results_file"])
                if not results_file.endswith(".csv"):
                    results_file = results_file + ".csv"
                result.to_csv(
                    cmd_pointer.workspace_path(cmd_pointer.settings["workspace"].upper()) + "/" + results_file,
                    index=False,
                )

        except Exception as e:
            print(e)
            result = response_result

        if isinstance(result, dict):
            if "error" in result:
                run_error = "Request Error:\n"
                for key, value in result["error"].items():
                    run_error = run_error + f"- <cmd>{key}</cmd> : {value}\n  "
                return output_text(run_error)

    except Exception as e:
        run_error = "HTTP Request Error:\n"

        return output_error(run_error + "\n" + str(e))

    return result


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

    auth_group = quoted_string | py.Word(py.alphanums + "_")
    service_name = quoted_string | py.Word(py.alphanums + "_")

    a_s = py.CaselessKeyword("as")
    describe = py.CaselessKeyword("describe")
    remote = py.CaselessKeyword("remote")
    auth = py.CaselessKeyword("auth")
    group = py.CaselessKeyword("group")
    _with = py.CaselessKeyword("with")
    add = py.CaselessKeyword("add")
    remove = py.CaselessKeyword("remove")
    to = py.CaselessKeyword("to")

    # catalog service
    using_keyword = py.CaselessKeyword("USING").suppress()
    quoted_identifier = py.QuotedString("'", escChar="\\", unquoteResults=True)
    parameter = py.Word(py.alphas, py.alphanums + "-_") | quoted_identifier
    value = py.Word(py.alphanums + "-_") | quoted_identifier
    param_value_pair = py.Group(parameter + py.Suppress("=") + value)
    using_clause = py.Optional(
        using_keyword + py.Suppress("(") + py.Optional(py.OneOrMore(param_value_pair))("params") + py.Suppress(")")
    )

    statements.append(py.Forward(model + auth + _list)("list_auth_services"))
    help.append(
        help_dict_create(
            name="model auth list",
            category="Model",
            command="model auth list",
            description="show authentication group mapping",
        )
    )

    statements.append(
        py.Forward(model + auth + add + group + auth_group("auth_group") + _with + quoted_string("api_key"))(
            "add_service_auth_group"
        )
    )
    help.append(
        help_dict_create(
            name="model auth add group",
            category="Model",
            command="model auth add group '<auth_group>'|<auth_group> with '<api_key>'",
            description="add an authentication group for model services to use",
        )
    )

    statements.append(py.Forward(model + auth + remove + group + auth_group("auth_group"))("remove_service_auth_group"))
    help.append(
        help_dict_create(
            name="model auth remove group",
            category="Model",
            command="model auth remove group '<auth_group>' | <auth_group>",
            description="remove an authentication group",
        )
    )

    statements.append(
        py.Forward(model + auth + add + service + service_name("service_name") + to + group + auth_group("auth_group"))(
            "attach_service_auth_group"
        )
    )
    help.append(
        help_dict_create(
            name="model auth add service",
            category="Model",
            command="model auth add service '<service_name>'|,service_name> to group '<auth_group>'|<auth_group>",
            description="Attach an authentication group to a model service",
        )
    )

    statements.append(
        py.Forward(model + auth + remove + service + service_name("service_name"))("detach_service_auth_group")
    )
    help.append(
        help_dict_create(
            name="model auth remove service",
            category="Model",
            command="model auth remove service '<service_name>'|<service_name>",
            description="Detatch an authentication group from a model service",
        )
    )

    statements.append(py.Forward(model + service + status)("model_service_status"))
    help.append(
        help_dict_create(
            name="model service status",
            category="Model",
            command="model service status",
            description="Get the status of currently cataloged services",
        )
    )

    statements.append(py.Forward(model + service + describe + (service_name)("service_name"))("model_service_config"))
    help.append(
        help_dict_create(
            name="model service describe",
            category="Model",
            command="model service describe '<service_name>'|<service_name>",
            description="get the configuration of a service",
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

    statements.append(py.Forward(uncatalog + model + service + service_name("service_name"))("uncatalog_model_service"))
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
            + using_clause
        )("catalog_add_model_service")
    )
    help.append(
        help_dict_create(
            name="catalog Model service",
            category="Model",
            command="catalog model service from (remote) '<path> or <github> or <service_url>' as  '<service_name>'|<service_name>   USING (<parameter>=<value> <parameter>=<value>)",
            description="""catalog a model service from a path or github or remotely from an existing OpenAD service.
(USING) optional headers parameters for communication with service backend.
If you are cataloging a service using a model defined in a directory, provide the absolute <cmd> <path> </cmd> of that directory in quotes.

The following options require the <cmd>remote</cmd> option be declared.

If you are cataloging a service using a model defined in github repository, provide the absolute <cmd> <github> </cmd> of that github directory quotes.

If you are cataloging a remote service on a ip address and port provide the remote services ipaddress and port in quoted string e.g. <cmd>'0.0.0.0:8080'</cmd>

<cmd>service_name</cmd>: this is the name of the service as you will define it for your usage. e.g <cmd>prop</cmd> short for properties. 

USING Parameters:

If using a hosted service the following parameters must be supplied:
-<cmd>Inference-Service</cmd>: this is the name of the inference service that is hosted, it is a required parameter if cataloging a remote service.
An authorization parameter is always required if cataloging a hosted service, either Auhtorisation group (<cmd>auth_group</cmd>) or Authorisation bearer_token/api_key (<cmd>Authorization</cmd>):
-<cmd>auth_group</cmd>: this is the name of an authorization group which contains the api_key linked to the service access. This can only be used if <cmd>Authorization</cmd> is not also defined.
OR
-<cmd>Authorization</cmd>: this parameter is designed to be used when a <cmd>auth_group</cmd> is not defined.

Example:

Skypilot Deployment
-<cmd>catalog model service from 'git@github.com:acceleratedscience/generation_inference_service.git' as 'gen'</cmd>

Service using a authentication group 
-<cmd>catalog model service from remote '<service_url>' as  molf  USING (Inference-Service=molformer  )</cmd>
<cmd> model auth add service 'molf' to group 'default'</cmd>

Single Authorisation Service
-<cmd>openad catalog model service from remote '<service_URL>' as 'gen' USING (Inference-Service=generation Authorization='<api_key>')</cmd>

Catalog a remote service shared with you:
-<cmd>catalog model service from remote 'http://54.235.3.243:30001' as gen</cmd>""",
        )
    )

    statements.append(
        py.Forward(
            model + service + up + service_name("service_name") + py.Optional(py.CaselessKeyword("NO_GPU")("no_gpu"))
        )("service_up")
    )
    help.append(
        help_dict_create(
            name="Model up",
            category="Model",
            command="model service up '<service_name>'|<service_name> [no_gpu]}",
            description="""launches a cataloged model service when it was cataloged as a self managed service from a directory or github repository.
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
            + service_name("service_name")
            + py.Optional(py.CaselessKeyword("NO_GPU")("no_gpu"))
        )("local_service_up")
    )
    help.append(
        help_dict_create(
            name="Model local up",
            category="Model",
            command="model service local up '<service_name>'|<service_name> ",
            description="""Launches a model service locally.

            Example:
              <cmd> model service local up gen</cmd>

             """,
        )
    )

    statements.append(py.Forward(model + service + down + service_name("service_name"))("service_down"))
    help.append(
        help_dict_create(
            name="Model down",
            category="Model",
            command="model service down '<service_name>'|<service_name>",
            description="Bring down a model service  \n Examples: \n\n<cmd>model service down gen</cmd> \n\n<cmd>model service down 'gen'</cmd> ",
        )
    )

    statements.append(
        py.Forward(
            py.CaselessKeyword("get")
            + model
            + py.CaselessKeyword("service")
            + service_name("service_name")
            + py.CaselessKeyword("result")
            + quoted_string("request_id")
        )("get_model_service_result")
    )
    help.append(
        help_dict_create(
            name="Get Model Service Result",
            category="Model",
            command="get model service '<service_name>'|<service_name> result '<result_id>' ",
            description="retrieves a result from a model service  \n Examples: \n\n<cmd>get model service myservier result 'wergergerg'  ",
        )
    )
