import os
import time
import shutil
import pytest
from pandas import DataFrame
from pandas.io.formats.style import Styler
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


@pytest.mark.parametrize(
    "from_path, service_name",
    [
        (
            "git@github.com:acceleratedscience/property_inference_service.git",
            random_name(),
        ),
        (
            "https://github.com/acceleratedscience/property_inference_service.git",
            random_name(),
        ),
        (f"/tmp/{random_name()}", random_name()),
    ],
)
def test_catalog_model_service(run_openad, from_path, service_name):
    """test run as user for adding a model service"""
    # setup
    cmd_add_service = f"catalog model service from '{from_path}' as '{service_name}'"
    cmd_remove_service = f"uncatalog model service '{service_name}'"
    if "/tmp" in from_path and not os.path.exists(from_path):
        # create tmp dir if testing local path
        os.makedirs(from_path)  # create tmp project to load into service
    # test
    result_add_service = run_openad(cmd_add_service)
    if "/tmp" in from_path and os.path.exists(from_path):
        # delete tmp dir if testing local path
        shutil.rmtree(from_path)  # immediately remove tmp project
    # validate service add
    with Dispatcher as test:
        assert result_add_service and service_name in test.list()
    # validate service remove
    result_remove_service = run_openad(cmd_remove_service)
    with Dispatcher as test:
        assert result_remove_service and service_name not in test.list()


# gloabal name to add and remove from test
up_service_name = random_name()


@pytest.mark.dependency()
def test_model_service_up(run_openad, dispatcher):
    """Test bringing up a model service"""
    # setup
    service_name = up_service_name
    from_path = "git@github.com:acceleratedscience/property_inference_service.git"
    cmd_model_status = "model service status"
    cmd_model_up = f"model service up '{service_name}'"
    cmd_add_service = f"catalog model service from '{from_path}' as '{service_name}'"
    # test and validate
    result_add_service = run_openad(cmd_add_service)
    assert result_add_service
    result_model_up = run_openad(cmd_model_up)
    # in case it fails dont want to leave hanging service
    if not result_model_up:
        with dispatcher() as cleanup:
            cleanup.remove_service(service_name)
    assert result_model_up
    time.sleep(1)
    timeout = 60 * 15  # 15 minutes
    start_t = time.time()
    status = "DOWN"
    while (time.time() - start_t) < timeout:
        if status == "PENDING":
            break
        df: DataFrame = run_openad(cmd_model_status)
        status = df[df["Service"].isin([service_name])]["Status"].iloc[0]
        print(status)
        time.sleep(10)
    assert status == "PENDING" or status == "READY"


@pytest.mark.dependency(depends=["test_model_service_up"])
def test_model_service_down(run_openad):
    """bring down service and uncatalog model service"""
    service_name = up_service_name
    cmd_model_down = f"model service down '{service_name}'"
    cmd_remove_service = f"uncatalog model service '{service_name}'"
    run_openad(cmd_model_down)
    run_openad(cmd_remove_service)
    with Dispatcher as test:
        assert service_name not in test.list()
