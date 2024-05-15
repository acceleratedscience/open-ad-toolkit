import json
from typing import Optional
from pydantic import BaseModel
from servicing import Dispatcher, UserProvidedConfig

### models are based off the servicing.bin serialized json ###


class ServiceData(BaseModel):
    remote_service: bool = False
    remote_endpoint: str = ""
    remote_up: bool = False


class ServiceConfig(BaseModel):
    """maps to UserProvidedConfig"""

    port: Optional[int] = None
    replicas: Optional[int] = None
    cloud: Optional[str] = None
    workdir: Optional[str] = None
    data: Optional[ServiceData] = ServiceData()
    disk_size: Optional[int] = None
    cpu: Optional[str] = None
    memory: Optional[str] = None
    accelerators: Optional[str] = None
    setup: Optional[str] = None
    run: Optional[str] = None

    def convert(self) -> UserProvidedConfig:
        """convert to UserProvidedConfig"""
        data = self.model_dump()
        data["data"] = json.dumps(self.data.model_dump())
        return UserProvidedConfig(**data)


if __name__ == "__main__":
    x = ServiceConfig()
    print(x.model_dump())
