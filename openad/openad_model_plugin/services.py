import datetime
from typing import Any, List, Optional, Callable, Tuple, Dict
from typing_extensions import Self
from servicing import Dispatcher, UserProvidedConfig
from openad.helpers.output import output_error, output_warning
import json
import os
from subprocess import run
import shlex
import requests
import time


class ModelServiceUniqueLocation(Dispatcher):
    """
    ModelService is a class that represents the servicing library
    """

    def __init__(self, cache: str | None = None) -> None:
        self.cache = cache
        self.name = "test"

        try:
            self.load(cache)
        except:
            self.create(cache)

    def save_on_exit(func: Callable[[], Any]) -> Any:
        """
        decorator to save the dispatcher services
        """

        def wrapper(self, *args, **kwargs) -> Any:
            print("> ---entering wrapper---")
            result = func(self, *args, **kwargs)
            print("> ---wrapper---", self)
            self.save()  # save the dispatcher state
            return result

        return wrapper

    def __call__(self, cache: Optional[str] = "") -> Self:
        if cache:
            self.cache = cache
        return self

    def __enter__(self) -> Self:
        self.load()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.save()
        self.cache = None

    def save(self, location: str | None = None) -> None:
        if location:
            # save to specified cache
            return super().save(location)
        elif self.cache:
            # save to saved cache
            return super().save(self.cache)
        else:
            # save to default cache
            return super().save()

    def create(self, location: str | None = None) -> None:
        """same signiture as save() but helps logic flow"""
        self.save(location)

    def load(self, location: str | None = None) -> None:
        model_cache = ""
        if location:
            # load specified cache
            model_cache = location
        elif self.cache:
            # load saved cache
            model_cache = self.cache
        else:
            # load default cache
            model_cache = None
        try:
            super().load(model_cache)
        except:
            # create and load services.bin file for namespace
            self.create(model_cache)
            self.load(model_cache)

    def add_service(self, name: str, config: UserProvidedConfig | None = None) -> None:
        return super().add_service(name, config)

    def remove_service(self, name: str) -> None:
        (f"removing service: {name}")
        return super().remove_service(name)

    def up(self, name: str) -> None:
        return super().up(name)

    def down(self, name: str) -> None:
        return super().down(name)

    def status(self, name: str, pretty: bool | None = None) -> Dict[str, Any]:
        return json.loads(super().status(name, pretty))

    def get_short_status(self, name: str):
        status = self.status(name)
        return {"up": status.get("up"), "url": status.get("url")}

    def get_url(self, name: str) -> str:
        return self.status(name).get("url")

    def get_services(self) -> List[Tuple[str, bool]]:
        """
        List the status of all the services"""
        details = []
        for service in self.list():
            status = self.status(service)
            details.append((service, status.get("up")))
        return details

    @property
    def servicer_apis(self) -> List[str]:
        """Gets available services during runtime

        Returns:
            List[str]: services
        """
        return list(filter(lambda x: not x.startswith("__"), dir(self)))


class ServiceFileLoadError(Exception):
    """Raises error if a service file could not load

    Args:
        Exception (_type_): _description_
    """

    def __init__(self, *args: object) -> None:
        super().__init__(*args)


class ModelService(Dispatcher):
    def __init__(self, *args, **kwargs) -> None:
        # search for previous running services
        self.load(location=kwargs.get("location"), update_status=kwargs.get("update_status"))
        super().__init__()

    def load(self, location: str | None = None, update_status: bool | None = False):
        """load a config. if it doesnt exist auto create it"""
        try:
            super().load(location=location, update_status=update_status)
            # output_success("loaded services")
        except:
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
                    output_warning(f"New services file created. Old services backed up to: {backup_location}")
                else:
                    # error
                    raise ServiceFileLoadError(
                        "Shutdown running services and try again or revert to an older version of openad"
                    )
            else:
                # ok. create file
                self.save(location=location)

    def __call__(self, *args: Any, **kwargs: Any) -> Self:
        # always does a load() but can optionally update servicer threads
        self.load(location=kwargs.get("location"), update_status=kwargs.get("update_status"))
        # TODO: load remote services here?
        return self

    def __enter__(self):
        # only does a load() if you call constructor method
        # like mymodelservice()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.save()

    def check_service_up(self, address: str, resource: str = "/health", timeout: float = 1.0) -> int | str:
        """ping the host address to see if service is up

        Args:
            address (str): _description_
            resource (str, optional): _description_. Defaults to "/health".
            timeout (float, optional): _description_. Defaults to 1.0.

        Returns:
            bool: _description_
        """
        start_time = time.time()
        # Ping the endpoint until timeout has elapsed
        up = False
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
            except requests.exceptions.RequestException as e:
                # Check if timeout has elapsed
                if time.time() - start_time >= timeout:
                    break
        return up

    def load_extra_data(self, name: str) -> Dict[str, Any]:
        """Returns data if data field in UserProvidedConfig is the only available field"""
        status = self.status(name).get("data")
        if status["data"]:
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
        return url

    def get_short_status(self, name: str) -> Dict[str, Any]:
        """Get service url and up status for local services and remote services

        Args:
            name (str): name of service

        Returns:
            Dict[str, Any]: {"up", "url"}
        """
        status = self.status(name)
        extra_data = self.load_extra_data(name)
        ret_status = {"is_remote": False}
        # alternative data exists
        if extra_data and extra_data.get("remote_endpoint"):
            url = extra_data.get("remote_endpoint")
            ret_status["url"] = extra_data.get("remote_endpoint")
            ret_status["up"] = self.check_service_up(url)
            ret_status["is_remote"] = True
        # use service data
        if status.get("url"):
            url = status.get("url")
        if status.get("up"):
            ret_status["up"] = bool(status.get("up"))
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

    def get_user_provided_config(self, name: str) -> UserProvidedConfig:
        status = self.status(name)
        return UserProvidedConfig(**status["data"])

    def get_config_as_dict(self, name: str) -> dict:
        return {**self.status(name)}

    def __get_build_step_count(self, name: str):
        dock_list = [
            "ADD",
            "ARG",
            "CMD",
            "COPY",
            "ENTRYPOINT",
            "ENV",
            "EXPOSE",
            "FROM",
            "HEALTHCHECK",
            "LABEL",
            "MAINTAINER",
            "ONBUILD",
            "RUN",
            "SHELL",
            "STOPSIGNAL",
            "USER",
            "VOLUME",
            "WORKDIR",
        ]
        status = self.status(name)
        workdir = status["template"]["workdir"]
        print(json.dumps(workdir, indent=2))
        with open(workdir + "/Dockerfile", "r") as f:
            dockerfile = [line.strip().split(" ", maxsplit=1)[0] for line in f.readlines()]
        # get a count for each dock in dockerfile
        total_dock_steps = 0
        for i in dockerfile:
            if i in dock_list:
                total_dock_steps += 1
        return total_dock_steps

    def get_build_log_completion(self, name: str):
        t_step = self.__get_build_step_count(name)
        cmd = shlex.split(f"sky serve logs {name} 1")
        print(cmd)
        print(run(cmd, capture_output=True))
        # run a subprocess to print output from stdout


class ServiceManager:
    def __init__(self) -> None:
        self.services = {}

    def __call__(self, location: str | None = None) -> ModelService:
        name = location
        if location is None:
            name = "default"
        if location in self.services:
            return self.services[name]
        else:
            service = ModelService(location)
            self.services.update({name: service})
            return self.services[name]


class DispatchManager:
    def __init__(self) -> None:
        self.head_dispatcher = Dispatcher()
        self.services = {}  # name: obj
        try:
            self.head_dispatcher.load()
        except:
            self.head_dispatcher.save()
        # load all services
        for name in self.head_dispatcher.list():
            self.catalog_service(name)

    def __call__(self, service_name: str | None = None) -> ModelService:
        print(f"dispatching {service_name}")
        return self.services.get(service_name)

    def catalog_service(self, service_name: str):
        self.services[service_name] = ModelService(service_name)

    def uncatalog_service(self, service_name: str):
        service = self.services.pop(service_name)
        del service


if __name__ == "__main__":
    dispatcher1 = Dispatcher()
    dispatcher1.load()
    print(dispatcher1.list())
    # print(json.dumps(dispatcher1.status(dispatcher1.list()[0]), indent=2))
