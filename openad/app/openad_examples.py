"""Populates Examples Directory"""
import os
import shutil
from openad.app.global_var_lib import _repo_dir


def openad_create_examps():
    """creates or overwrites example directory"""

    try:
        try:
            shutil.rmtree(os.path.expanduser("~/openad_notebooks"), ignore_errors=True)
        except:  # nothing to delete
            pass

        shutil.copytree(
            _repo_dir + "/../notebooks",
            os.path.abspath(os.path.expanduser("~/openad_notebooks")),
            dirs_exist_ok=True,
            symlinks=True,
        )

    except Exception as e:
        print("Error: Unable to copy samples")
        print(e)
