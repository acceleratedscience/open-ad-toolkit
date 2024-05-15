import datetime
from typing import Any, Dict
from typing_extensions import Self
from servicing import Dispatcher, UserProvidedConfig
from openad.helpers.output import output_error, output_warning
import json
import os
import requests
import time
from openad.openad_model_plugin.models import ServiceConfig
from openad.openad_model_plugin.utils import get_logger


logger = get_logger(__name__)


class ServiceFileLoadError(Exception):
    """Raises error if a service file could not load"""

    def __init__(self, *args: object) -> None:
        super().__init__(*args)


class Service:
    """class to simplify service data fetching"""

    def __init__(self, name, location) -> None:
        self.name: str = name
        self.location = location
        self.dispatcher: Dispatcher = Dispatcher(skip_sky_validation=True)
        self.dispatcher.load()
        self.status_copy: dict = self.status

    @property
    def status(self) -> dict:
        # return old copy if service is removed but need its data
        try:
            logger.debug(f"fetching status '{self.name=}' | {self.location=}")
            # TODO: maybe cache status per context manager?
            return json.loads(self.dispatcher.status(self.name))
            return self.dispatcher.status(self.name)
        except Exception as e:
            logger.error(f"'{self.name=}' | {str(e)}. returning cached status")
            return self.status_copy

    @property
    def userprovidedconfig(self) -> UserProvidedConfig:
        return UserProvidedConfig(**self.status.get("data"))

    @property
    def serviceconfig(self) -> ServiceConfig:
        data = self.config
        data["data"] = json.dumps(self.config_data)
        return ServiceConfig(**data)

    @property
    def config(self) -> dict:
        return self.status.get("data")

    @property
    def config_data(self) -> dict:
        data = self.config.get("data", {}) or json.dumps({})
        return json.loads(data)

    @property
    def is_remote(self) -> bool:
        return self.config_data.get("remote_service", False)

    @property
    def is_alive(self) -> bool:
        if self.is_remote:
            return self.check_service_up()
        return self.status.get("up", False)
        # return self.status.get("up", False) or self.check_service_up()

    @property
    def template(self) -> dict:
        return self.status.get("template", {})

    def check_service_up(self, resource: str = "/health", timeout: float = 1.0) -> bool:
        """ping the host address to see if service is up

        Args:
            resource (str, optional): _description_. Defaults to "/health".
            timeout (float, optional): _description_. Defaults to 1.0.

        Returns:
            bool: _description_
        """
        start_time = time.time()
        # Ping the endpoint until timeout has elapsed
        up = False
        address = self.url  # reduces calls to function
        while True:
            try:
                # Make a GET request to the endpoint
                response = requests.get(address + resource, timeout=0.2)
                if response.status_code == 200:
                    up = True
                    break
                # Check if timeout has elapsed
                if time.time() - start_time >= timeout:
                    break
            except requests.exceptions.RequestException:
                # Check if timeout has elapsed
                if time.time() - start_time >= timeout:
                    break
        return up

    @property
    def url(self) -> str:
        url = ""
        # check if remote url
        if self.is_remote:
            url = self.config_data.get("remote_endpoint")
        elif self.status.get("url"):
            url = self.status.get("url")
        # add url prefix
        if url and "http" not in url:
            url = "http://" + url
        return url


class ModelService(Dispatcher):
    def __init__(
        self,
        location: str | None = None,
        update_status: bool | None = False,
        skip_sky_validation: bool = False,
    ) -> None:
        self.services: Dict[str, Service] = dict()
        # search for previous running services
        self.load(location=location, update_status=update_status)
        super().__init__()

    def save(self, location: str | None = None) -> None:
        # check if service deleted by another window
        # dispatch = Dispatcher(skip_sky_validation=True)
        # if sorted(dispatch.list()) != self.list():
        #     logger.debug("services changed!. not saving")
        #     return
        logger.debug(f"saved config {location=}")
        return super().save(location)

    def load(
        self,
        location: str | None = None,
        update_status: bool | None = False,
        extend: bool = True,
    ):
        """load a config. if it doesnt exist auto create it"""
        try:
            # load services
            super().load(extend=extend, location=location, update_status=update_status)
            logger.debug(f"loaded config {location=} {update_status=} {extend=}")
            # create service objects
            for name in self.list():
                if name not in self.services:
                    self.update_service_map(name, location=location)
            logger.debug(f"config objects  | num={len(self.list())} | {self.list()}")
            logger.debug(
                f"Service objects | num={len(self.services)} | {list(self.services.keys())}"
            )

        except Exception as e:
            output_error(str(e))
            # use default directory
            if location is None:
                location = os.path.expanduser("~/.servicing")
            # check if services file exists in path
            if os.path.exists(location + "/services.bin"):
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
                    output_warning(
                        f"New services file created. Old services backed up to: {backup_location}"
                    )
                else:
                    # error
                    raise ServiceFileLoadError(
                        "Shutdown running services and try again or revert to an older version of openad"
                    )
            else:
                # ok. create file
                self.save(location=location)

    def __call__(
        self, location: str | None = None, update_status: bool = False
    ) -> Self:
        # always does a load() but can optionally update servicer threads
        self.load(location=location, update_status=update_status)
        # TODO: load remote services here?
        return self

    def __enter__(self):
        # only does a load() if you call constructor method
        # with mymodelservice() as service:
        logger.debug("---entering context---")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.save()
        logger.debug("---exiting context---")

    def up(
        self, name: str, skip_prompt: bool | None = None, gpu_disable: bool = False
    ) -> None:
        service = self.services.get(name)
        logger.debug(f"starting service '{name}'")
        # TODO: update openad.cfg file for resource state
        if gpu_disable:
            logger.debug("gpu flag found")
            # remove gpu from config
            if service.config.get("accelerators"):
                try:
                    logger.debug("creating temporary service without gpu")
                    # get resource values for service
                    tmp_data = service.serviceconfig
                    tmp_data.accelerators = None
                    gpu_disable_config = ServiceConfig(**tmp_data)
                    # replace current config
                    self.remove_service(name)
                    self.add_service(name, config=gpu_disable_config)
                    # start service without gpu
                    super().up(name, skip_prompt)
                except Exception as e:
                    # put back old config
                    self.remove_service(name)
                    self.add_service(name, config=service.userprovidedconfig)
                    raise e
                finally:
                    # put back old config
                    self.remove_service(name)
                    self.add_service(name, config=service.userprovidedconfig)
                    logger.debug(f"reverted '{name}' to original configuration")
                    return
            else:
                output_warning("service already not configured for gpu")
        return super().up(name, skip_prompt)

    def down(
        self, name: str, skip_prompt: bool | None = None, force: bool | None = None
    ) -> None:
        logger.debug(f"shutting down service '{name}' | {skip_prompt=} {force=}")
        service = self.services.get(name)
        if service.is_remote:
            logger.debug(f"service '{name}' is remote service. exiting.")
            return
        # stop service
        super().down(name, skip_prompt, force)
        # reinitialize service
        self.reinit_service(name, service.userprovidedconfig)

    def reinit_service(
        self,
        name: str,
        config: ServiceConfig | UserProvidedConfig | None = None,
        remote_endpoint: str | None = None,
    ):
        """same signiture as add_service but removes and readds service"""
        logger.debug(f"reinit '{name}'")
        if name in self.list():
            # remove the service without removing self.services[name] object
            super().remove_service(name)
            self.add_service(name, config=config)
            logger.debug(f"'{name}' reinitialized.")

    def remove_service(self, name: str) -> None:
        if name in self.list():
            logger.debug(f"removing service '{name}'")
            self.services.pop(name)  # remove service from map
            super().remove_service(name)

    def add_service(
        self,
        name: str,
        location: str,
        config: ServiceConfig | UserProvidedConfig | None = None,
        remote_endpoint: str | None = None,
    ) -> None:
        logger.debug(
            f"adding service '{name}' | type={type(config)} {remote_endpoint=}"
        )
        if isinstance(config, UserProvidedConfig):
            super().add_service(name, config)
        elif not config:
            config = ServiceConfig()
            if remote_endpoint:
                config.data.remote_service = True
                config.data.remote_endpoint = remote_endpoint
        else:
            super().add_service(name, config.convert())
        # have to add after service added
        self.update_service_map(name, location, config)

    def update_service_map(self, name: str, location: str):
        if name not in self.services:
            self.services[name] = Service(name, location)
            logger.debug(f"creating object Service({name})")

    def get_url(self, name: str) -> str:
        return self.services.get(name).url

    def get_short_status(self, name: str) -> Dict[str, Any]:
        """Get service url and up status for local services and remote services

        Args:
            name (str): name of service

        Returns:
            Dict[str, Any]: {"up", "url"}
        """
        service = self.services.get(name)
        ret_status = {
            "is_remote": service.is_remote,
            "up": service.is_alive,
            "url": service.url,
        }
        return ret_status

    def status(self, name: str, pretty: bool | None = None) -> Dict[str, Any]:
        """Loads status as json object

        Args:
            name (str): service name
            pretty (bool | None, optional): format indent. Defaults to None.

        Returns:
            Dict[str, Any]: data
        """
        return json.loads(super().status(name, pretty))

    def get_config_as_dict(self, name: str) -> dict:
        return self.services.get(name).status.copy()

    def get(self, name: str) -> Service | None:
        return self.services.get(name)


if __name__ == "__main__":
    dispatcher1 = ModelService(skip_sky_validation=True)
    dispatcher1.load()
    print(dispatcher1.list())
    print(dispatcher1.services["remote"].is_remote)
    print(json.dumps(dispatcher1.status(dispatcher1.list()[0]), indent=2))
