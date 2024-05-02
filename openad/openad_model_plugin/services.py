from typing import Any, List, Optional, Callable, Tuple, Dict
from typing_extensions import Self
from servicing import Dispatcher, UserProvidedConfig
from openad.helpers.output import output_error, output_text, output_success
import json
import os
from subprocess import run
import shlex


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


class ModelService(Dispatcher):
    def __init__(self, location: str = None, update_status: bool = True) -> None:
        super().__init__()
        # search for previous running services
        self.load(location=location, update_status=update_status)

    def load(self, location: str = None, update_status: bool = False):
        """load a config. if it doesnt exist auto create it"""
        try:
            super().load(location=location, update_status=update_status)
        except:
            self.save(location=location)

    def __enter__(self, name: str = None):
        self.name = name
        self.load()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.save()

    def get_short_status(self, name: str):
        status = self.status(name)
        return {"up": bool(status.get("up")), "url": status.get("url")}

    def status(self, name: str, pretty: bool | None = None) -> Dict[str, Any]:
        return json.loads(super().status(name, pretty))

    def get_config(self, name: str):
        status = self.status(name)
        return UserProvidedConfig(**status["data"])

    def down(self, name: str, force: bool = None):
        super().down(name, force)

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
    print(json.dumps(dispatcher1.status(dispatcher1.list()[0]), indent=2))
