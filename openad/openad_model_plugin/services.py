from typing import Any, List, Optional, Callable, Tuple, Dict
from typing_extensions import Self
from servicing import Dispatcher, UserProvidedConfig
import json


class ModelService(Dispatcher):
    """
    ModelService is a class that represents the servicing library
    """

    def __init__(self) -> None:
        self.cache = ""
    
    def __call__(self, cache: Optional[str]="") -> Self:
        if cache:
            self.cache = cache
        return self
    
    def __enter__(self) -> Self:
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
            self.cache = ""  # reset cache dir
        else:
            print(f"saving from default")
            self.save()

    def set_cache(self, path: str) -> None:
        self.cache = path
    
    def download_model(self, name:str, url: str, model_dir: str):
        pass
    
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

    def load(self, location: str | None = None) -> None:
        if location:
            # load specified cache
            return super().load(location)
        elif self.cache:
            # load saved cache
            return super().load(self.cache)
        else:
            # load default cache
            return super().load()

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
    def add_service(self, name: str, config: UserProvidedConfig | None = None) -> None:
        return super().add_service(name, config)

    @save_on_exit
    def remove_service(self, name: str) -> None:
        return super().remove_service(name)

    @save_on_exit
    def up(self, name: str) -> None:
        return super().up(name)

    @save_on_exit
    def down(self, name: str) -> None:
        return super().down(name)

    def status(self, name: str, pretty: bool | None = None) -> Dict[str, Any]:
        return json.loads(super().status(name, pretty))

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


if __name__ == "__main__":
    t = ModelService()
    print(t.servicer_apis)
