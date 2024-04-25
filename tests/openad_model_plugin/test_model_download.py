from openad.openad_model_plugin.catalog_model_services import retrieve_model
import pytest
import os
import shutil
import random


@pytest.mark.parametrize(
    'from_path',
    [
        "git@github.com:acceleratedscience/property_inference_service.git",
        "git@github.com:acceleratedscience/generation_inference_service.git",
        "https://github.com/acceleratedscience/property_inference_service.git",
        pytest.param("/tmp/testingopenad", marks=pytest.mark.xfail(reason="path doesnt exist")),
    ]
)
def test_download_with_git(from_path):
    # setup
    local_to_path = "/tmp/tests"
    # test
    is_success, err = retrieve_model(from_path, local_to_path)
    assert is_success, str(err)
    shutil.rmtree(os.path.join(local_to_path))


def test_invalid_local_path():
    # setup
    local_from_path = "/tmp/tests" + bin(random.getrandbits(10))
    local_to_path = "/tmp/tests" + bin(random.getrandbits(10))
    is_success, err = retrieve_model(local_from_path, local_to_path)
    assert not is_success
