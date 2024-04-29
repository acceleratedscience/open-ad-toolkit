import pytest
from openad.openad_model_plugin.services import ModelService
from openad.openad_model_plugin.catalog_model_services import *

@pytest.fixture(scope="session")
def model_service():
    return ModelService()

def test_service_init(model_service):
    assert isinstance(model_service.list(), list)

@pytest.mark.skip
def test_model_service_status():
    cmd_pointer, parser = "", ""
    # print(model_service_status(cmd_pointer, parser))