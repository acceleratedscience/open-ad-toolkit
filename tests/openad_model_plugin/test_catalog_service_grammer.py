import os
import time
import shutil
import pytest
from pandas import DataFrame
from openad.openad_model_plugin.catalog_model_services import Dispatcher
from tests.helpers import random_name


def test_model_service_status(run_openad):
    """get model status"""
    # setup
    cmd = "model service status"
    # test
    result = run_openad(cmd)
    # validate
    assert isinstance(result, DataFrame)


def test_catalog_model_service(run_openad):
    """test run as user for adding a model service"""
    # setup
    service_name = random_name()
    path = "/tmp/" + service_name
    cmd_add_service = f"catalog model service from '{path}' as '{service_name}'"
    cmd_remove_service = f"uncatalog model service '{service_name}'"
    if not os.path.exists(path):
        os.makedirs(path)  # create tmp project to laod into service
    # test
    result = run_openad(cmd_add_service)
    if os.path.exists(path):
        shutil.rmtree(path)  # immediately remove tmp project
    # validate service add
    with Dispatcher as test:
        assert service_name in test.list()
    # validate service remove
    result = run_openad(cmd_remove_service)
    with Dispatcher as test:
        assert service_name not in test.list()


@pytest.mark.skip
@pytest.mark.dependency()
def test_model_service_up(run_openad):
    """Test bringing up a model service"""
    # setup
    service_name = "gt4sd_prop"
    from_path = "git@github.com:acceleratedscience/property_inference_service.git"
    cmd_model_status = "model service status"
    cmd_model_up = f"model service up '{service_name}'"
    cmd_add_service = f"catalog model service from '{from_path}' as '{service_name}'"
    # test
    result1 = run_openad(cmd_add_service)
    result2 = run_openad(cmd_model_up)
    time.sleep(1)
    # validate
    timeout = 60 * 15  # 15 minutes
    start_t = time.time()
    status = "DOWN"
    while (time.time() - start_t) < timeout:
        if status == "PENDING":
            break
        df: DataFrame = run_openad(cmd_model_status)
        status = df[df['Service'].isin([service_name])]['Status'].iloc[0]
        print(status)
        time.sleep(10)
    assert status == "PENDING"


@pytest.mark.dependency(depends=['test_model_service_up'])
def test_model_service_down(run_openad):
    """bring down service and uncatalog model service"""
    service_name = "gt4sd_prop"
    cmd_model_down = f"model service down '{service_name}'"
    cmd_remove_service = f"uncatalog model service '{service_name}'"
    run_openad(cmd_model_down)
    run_openad(cmd_remove_service)
    with Dispatcher as test:
        assert service_name not in test.list()