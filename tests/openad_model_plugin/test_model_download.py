from openad.openad_model_plugin.catalog_model_services import retrieve_model
import pytest
import os
import shutil
import random


@pytest.fixture
def run_file_cleanup(tmpdir):
    """Fixture to execute asserts before and after a test is run"""
    # Setup: fill with any logic you want

    yield # this is where the testing happens

    # Teardown : fill with any logic you want
    if os.path.exists(tmpdir):
        print("removing dir:", tmpdir)
        shutil.rmtree(tmpdir)


@pytest.mark.parametrize(
    'from_path',
    [
        "git@github.com:acceleratedscience/servicing.git",
        "https://github.com/acceleratedscience/servicing.git",
        "/tmp/testing",
    ]
)
def test_download(from_path):
    # setup
    local_to_path = "/tmp/tests"
    # test
    is_success, err = retrieve_model(from_path, local_to_path)
    assert is_success, str(err)
    shutil.rmtree(os.path.join(local_to_path))


def test_invalid_path():
    # setup
    local_from_path = "/tmp/tests" + bin(random.getrandbits(10))
    local_to_path = "/tmp/tests" + bin(random.getrandbits(10))
    is_success, err = retrieve_model(local_from_path, local_to_path)
    assert not is_success