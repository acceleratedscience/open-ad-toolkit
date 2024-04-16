from typing import Any, List, Optional, Callable, Tuple, Dict
from servicing import Dispatcher, UserProvidedConfig

import json


class ModelService:
    """
    ModelService is a class that represents the servicing library
    """

    def __init__(self, cache: Optional[str]="") -> None:
        self.dispatcher = Dispatcher()
        self.cache = cache
    
    def __enter__(self):
        if self.cache:
            print(f"loading from cache: {self.cache}")
            self.load(self.cache)
        else:
            print(f"loading from default")
            self.load()
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.cache:
            print(f"saving from cache: {self.cache}")
            self.save(self.cache)
        else:
            print(f"saving from default")
            self.save()

    def set_cache(self, path: str) -> None:
        self.cache = path
    
    def download_model(self, name:str, url: str, model_dir: str):
        pass

    def save(self, location: Optional[str] = None):
        if location:
            # save to specified cache
            self.dispatcher.save(location)
        elif self.cache:
            # save to saved cache
            self.dispatcher.save(self.cache)
        else:
            # save to default cache
            self.dispatcher.save()

    def load(self, location: Optional[str] = None):
        if location:
            # load specified cache
            self.dispatcher.load(location)
        elif self.cache:
            # load saved cache
            self.dispatcher.load(self.cache)
        else:
            # load default cache
            self.dispatcher.load()

    def save_on_exit(func: Callable[[], Any]) -> Any:
        """
        decorator to save the dispatcher services
        """
        def wrapper(self, *args, **kwargs) -> Any:
            result = func(self, *args, **kwargs)
            self.save()  # save the dispatcher state
            return result
        return wrapper

    @save_on_exit
    def add_service(self, name: str,
                    config: Optional[UserProvidedConfig] = None) -> None:
        self.dispatcher.add_service(name, config)

    @save_on_exit
    def remove_service(self, name: str) -> None:
        self.dispatcher.remove_service(name)

    @save_on_exit
    def up(self, name: str) -> None:
        self.dispatcher.up(name)

    @save_on_exit
    def down(self, name: str) -> None:
        self.dispatcher.down(name)
    
    def list(self) -> List[str]:
        return self.dispatcher.list()
    
    def status(self, name: str, pretty: Optional[bool] = None) -> Dict[str,str]:
        return json.loads(self.dispatcher.status(name, pretty))

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
        return list(filter(lambda x: not x.startswith("__"), dir(self.dispatcher)))


if __name__ == "__main__":
    t = ModelService()
    print(t.servicer_apis)
