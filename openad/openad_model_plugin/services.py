import os
import time
import datetime
import json
import requests
from typing import Any, Dict
from typing_extensions import Self
from openad.helpers.output import output_error, output_warning
from servicing import Dispatcher, UserProvidedConfig
from openad.openad_model_plugin.utils import get_logger, LruCache
from functools import cache, lru_cache


logger = get_logger(__name__)


remote_services_cache = LruCache[dict](capacity=16)


class ServiceFileLoadError(Exception):
    """Raises error if a service file could not load

    Args:
        Exception (_type_): _description_
    """

    def __init__(self, *args: object) -> None:
        super().__init__(*args)


class ServiceFetchError(Exception):
    """Raises error when failed to fetch remote service definitions

    Raises:
        ServiceFileLoadError:
        e: failed to fetch remote service definitions
    """

    def __init__(self, *args: object) -> None:
        super().__init__(*args)


class ModelService(Dispatcher):
    def __init__(
        self,
        location: str = None,
        update_status: bool = False,
        skip_sky_validation: bool = True,
    ) -> None:
        logger.debug(f"dispatcher init | {location=} {update_status=} {skip_sky_validation=}")
        self.default_location = location
        self.load(location=location, update_status=update_status)
        super().__init__()

    def save(self, location: str | None = None) -> None:
        location = location or self.default_location
        logger.debug(f">> saving config << | {location=}")
        return super().save(location)

    def load(self, location: str | None = None, update_status: bool | None = False):
        """load a config. if it doesnt exist auto create it"""
        location = location or self.default_location
        logger.debug(f">> loading config << | {location=} {update_status=}")
        try:
            super().load(location=location, update_status=update_status)
            # output_success("loaded services")
        except:
            # use default directory
            if location is None:
                location = os.path.expanduser("~/.servicing")
            # check if services file exists in path
            if os.path.exists(os.path.join(location, "services.bin")):
                logger.warn(f"config already exists but could not load | {location=} {update_status=}")
                output_error("Error: unable to load services config")
                user_choice = input(
                    "Do you want to overwrite the current services?\nUse with caution this can lead to hanging services! (y/n) : "
                )
                if user_choice.strip().lower() == "y":
                    # rename directory to save as backup
                    backup_location = f"{location}.backup-{datetime.datetime.now().strftime('%Y-%m-%d-%H:%M:%S')}"
                    os.rename(location, backup_location)
                    # create new services directory
                    self.save(location=location)
                    output_warning(f"New services file created. Old services backed up to: {backup_location}")
                else:
                    # error
                    logger.error(f"creating config | {location=} {update_status=}")
                    raise ServiceFileLoadError(
                        "Shutdown running services and try again or revert to an older version of openad"
                    )
            else:
                # ok. create file
                logger.debug(f"creating config | {location=} {update_status=}")
                self.save(location=location)

    def __call__(self, location: str | None = None, update_status: bool = False) -> Self:
        if location or update_status:
            logger.debug(f"update load method | {location=} {update_status=}")
        # always does a load() but can optionally update servicer threads
        self.load(location=location, update_status=update_status)
        return self

    def __enter__(self):
        # only does a load() if you call constructor method
        # like mymodelservice()
        # logger.debug("--entering context--")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        # logger.debug("--exiting context--")
        pass

    def load_as_b64(self, b64: str) -> None:
        # No Implementation
        logger.debug(f"loading config from b64 string | {b64=}")
        super().load_as_b64(b64)
        self.save()

    def down(self, name: str, skip_prompt: bool | None = None, force: bool | None = None) -> None:
        logger.debug(f"stopping running service | {name=} {skip_prompt=} {force=}")
        super().down(name, skip_prompt, force)
        self.save()

    def remove_service(self, name: str) -> None:
        logger.debug(f"removing service | {name=}")
        super().remove_service(name)
        self.save()

    def add_service(self, name: str, config: UserProvidedConfig | None = None) -> None:
        logger.debug(f"adding service | {name=}")
        super().add_service(name, config)
        self.save()

    def up(self, name: str, skip_prompt: bool | None = None, gpu_disable: bool = False) -> None:
        # TODO: update openad.cfg file for resource state
        logger.debug(f"starting service | {name=} {skip_prompt=} {gpu_disable=}")
        if gpu_disable:
            service_config_dict = self.get_config_as_dict(name)
            old_config = self.get_user_provided_config(name)
            # remove gpu from config
            if service_config_dict.get("data")["accelerators"]:
                try:
                    # get resource values for service
                    tmp_data = {**service_config_dict["data"]}
                    tmp_data["accelerators"] = None
                    gpu_disable_config = UserProvidedConfig(**tmp_data)
                    # replace current config
                    self.remove_service(name)
                    self.add_service(name, config=gpu_disable_config)
                    # start service without gpu
                    super().up(name, skip_prompt)
                except Exception as e:
                    # put back old config
                    self.remove_service(name)
                    self.add_service(name, config=old_config)
                    raise e
                finally:
                    # put back old config
                    self.remove_service(name)
                    self.add_service(name, config=old_config)
                    return
            else:
                output_warning("service already not configured for gpu")
        super().up(name, skip_prompt)
        self.save()

    def check_service_up(self, address: str, resource: str = "/health", timeout: float = 2.0) -> int | str:
        """ping the host address to see if service is up

        Args:
            address (str): _description_
            resource (str, optional): _description_. Defaults to "/health".
            timeout (float, optional): _description_. Defaults to 1.0.

        Returns:
            bool: _description_
        """
        logger.debug(f"pinging service | {address=} {resource=} {timeout=}")
        up = False
        # Ping the endpoint until timeout has elapsed
        try:
            # Make a GET request to the endpoint
            response = requests.get(address + resource, timeout=timeout)
            if response.status_code == 200:
                up = True
            else:
                # got some error status code
                logger.debug(f"got unexpected status code of {response.status_code}")
        except requests.exceptions.RequestException:
            logger.debug(f"failed to reach {address=}")
        logger.debug(f"service endpoint status | {address+resource} {up=}")
        return up

    def load_extra_data(self, name: str) -> Dict[str, Any]:
        """Returns data if data field in UserProvidedConfig is the only available field"""
        status = self.status(name).get("data")
        if status["data"]:
            logger.debug(f"service has additional data | {name=}")
            return json.loads(status["data"])
        return dict()

    def get_url(self, name: str):
        status = self.status(name)
        extra_data = self.load_extra_data(name)
        url = ""
        # check if remote url
        if extra_data and extra_data.get("remote_endpoint"):
            url = extra_data.get("remote_endpoint")
        elif status.get("url"):
            url = status.get("url")
        # add url prefix
        if url and "http" not in url:
            url = "http://" + url
        logger.debug(f"formatted service url | {name=} {url=}")
        return url

    def get_short_status(self, name: str) -> Dict[str, Any]:
        """Get service url and up status for local services and remote services

        Args:
            name (str): name of service

        Returns:
            Dict[str, Any]: {"up", "url"}
        """
        # TODO: refactor
        status = self.status(name)
        extra_data = self.load_extra_data(name)
        # init default status data
        ret_status = {"is_remote": False, "url": "", "up": False}
        # alternative data exists
        if extra_data and extra_data.get("remote_endpoint"):
            # use remote service data
            url = extra_data.get("remote_endpoint")
            ret_status["url"] = extra_data.get("remote_endpoint")
            ret_status["up"] = self.check_service_up(url)
            ret_status["is_remote"] = True
        # use local service data
        if status.get("url"):
            ret_status["url"] = self.get_url(name)
        if status.get("up"):
            ret_status["up"] = bool(status.get("up"))
        elif not ret_status.get("is_remote") and ret_status.get("url"):
            # TODO: this should be fixed in servicing library
            # recheck if service is actually down
            ret_status["up"] = self.check_service_up(ret_status["url"])
        logger.debug(f"service info | {name=} {ret_status=}")
        return ret_status

    def get_remote_service_definitions(self, name: str) -> list | None:
        """retrieve remote service definitions. caches first result"""
        if remote_services_cache.get(name):
            # cache hit
            logger.debug(f"getting '{name}' from cache")
            return remote_services_cache.get(name)
        service_definitions = []
        service_data = self.get_short_status(name)
        if service_data.get("is_remote"):
            logger.debug(f"fetching remote service defs | {name=}")
            try:
                url = service_data.get("url")
                response = requests.get(url + "/service", timeout=10)
                if response.status_code == 200:
                    service_definitions = response.json()
            except requests.exceptions.RequestException as e:
                # could not get service defs. service not up or wrong url
                logger.error(f"could not fetch service defs | {name=} {str(e)}")
        if service_definitions:
            # insert when not None
            logger.debug(f"inserting '{name}' into cache")
            remote_services_cache.insert(name, service_definitions)
        return service_definitions

    def get_service_cache(self) -> LruCache[dict]:
        return remote_services_cache

    def status(self, name: str, pretty: bool | None = None) -> Dict[str, Any]:
        """Loads status as json object

        Args:
            name (str): service name
            pretty (bool | None, optional): format indent. Defaults to None.

        Returns:
            Dict[str, Any]: data
        """
        return json.loads(super().status(name, pretty))

    def get_user_provided_config(self, name: str) -> UserProvidedConfig:
        logger.debug(f"getting provider config | {name=}")
        status = self.status(name)
        return UserProvidedConfig(**status["data"])

    def get_config_as_dict(self, name: str) -> dict:
        logger.debug(f"getting config as dict | {name=}")
        return {**self.status(name)}


if __name__ == "__main__":
    dispatcher1 = ModelService()
    dispatcher1.load()
    print(dispatcher1.list())
    print(json.dumps(dispatcher1.status(dispatcher1.list()[0]), indent=2))
