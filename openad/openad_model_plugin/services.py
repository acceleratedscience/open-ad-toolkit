from typing import Any, List, Optional, Callable, Tuple
from servicing import Dispatcher, UserProvidedConfig
import json


class ModelService(Dispatcher):
    """
    ModelService is a class that represents the servicing library
    """

    def __init__(self, location: Optional[str|None]=None) -> None:
        # load the services on init
        if location:
            self.load(location)

    def save_on_exit(func: Callable[[], Any]) -> Any:
        """
        decorator to save the dispatcher services
        """
        def wrapper(self, *args, **kwargs) -> Any:
            result = func(self, *args, **kwargs)
            self.save()  # save the dispatcher state
            return result
        return wrapper

    def load(self, location: Optional[str] = None):
        super().load(location)
        return self  # chain methods

    @save_on_exit
    def add_service(self, name: str,
                    config: Optional[UserProvidedConfig] = None) -> None:
        super().add_service(name, config)

    @save_on_exit
    def remove_service(self, name: str) -> None:
        super().remove_service(name)

    @save_on_exit
    def up(self, name: str) -> None:
        super().up(name)

    @save_on_exit
    def down(self, name: str) -> None:
        super().down(name)

    def get_services(self) -> List[Tuple[str, bool]]:
        """
        List the status of all the services"""
        details = []
        for service in self.list():
            status = json.loads(self.status(service))
            details.append((service, status.get("up") ))
        return details

    @property
    def servicer_apis(self) -> List[str]:
        """Gets available services during runtime

        Returns:
            List[str]: services
        """
        return list(filter(lambda x: not x.startswith("__"), dir(self)))
