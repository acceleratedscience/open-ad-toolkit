import pytest
from openad.openad_model_plugin.services import ModelService

@pytest.fixture(scope="session")
def model_service():
    return ModelService()

def test_service_init(model_service):
    assert isinstance(model_service.list(), list)