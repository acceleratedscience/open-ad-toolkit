from typing import Any, List, Optional, Callable, Tuple, Dict
from typing_extensions import Self
from servicing import Dispatcher, UserProvidedConfig
import json


class ModelService(Dispatcher):
    """
    ModelService is a class that represents the servicing library
    """
    def __init__(self, cache: str | None = None) -> None:
        self.cache = cache
        self.name = "test"
    
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
    
    def __call__(self, cache: Optional[str]="") -> Self:
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
        print(f'removing service: {name}')
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
            details.append((service, status.get("up") ))
        return details

    @property
    def servicer_apis(self) -> List[str]:
        """Gets available services during runtime

        Returns:
            List[str]: services
        """
        return list(filter(lambda x: not x.startswith("__"), dir(self)))


class ServiceManager:
    def __init__(self) -> None:
        self.services = {}
    def __call__(self, location: str | None = None) -> ModelService:
        name = location
        if location is None: name = "default"
        if location in self.services:
            return self.services[name]
        else:
            service = ModelService(location)
            self.services.update({name: service})
            return self.services[name]


if __name__ == "__main__":
    import os
    location = os.path.expanduser("~/.openad_model_services/gt4sd_prop")
    servicer = ServiceManager()
    
    print(servicer().cache)
    with servicer(location) as model:
        print(model.cache)
    with servicer(location) as model:
        print(model.list())
