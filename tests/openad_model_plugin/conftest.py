import pytest
from openad.app.global_var_lib import GLOBAL_SETTINGS
from openad.app.main import RUNCMD
from openad.openad_model_plugin.services import ModelService


### global pytest fixtures ###


@pytest.fixture(scope="session")
def run_openad():
    """Run openad commands as user. mimics cli"""
    # set to return output
    GLOBAL_SETTINGS["display"] = "notebook"
    # return function call
    return RUNCMD().default


@pytest.fixture(scope="session")
def dispatcher():
    """
    Create a new model service instance that can be refrenced.
    Recommend using the running object from catalog_model_services.py (Dispatcher)
    """
    return ModelService()
