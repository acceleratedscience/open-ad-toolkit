import os
from jupyter_core import paths

try:
    import jupyterlab_server.config as jls
except:
    pass


def get_jupyter_lab_config():
    """pulls the jupyter lab settings as dictionary"""
    # Get the path to the Jupyter server configuration directory
    config_dir = ""
    try:
        config_dir = jls.jupyter_config_dir()
    except:
        try:
            config_dir = paths.jupyter_config_dir()
        except:
            config_dir = ""

    # Define the configuration file path

    config_file = os.path.join(config_dir, "jupyter_lab_config.py")

    jupyter_lab_settings = {}
    # Check if the configuration file exists
    if os.path.exists(config_file):
        with open(config_file, "r") as f:
            config_lines = f.readlines()
        for line in config_lines:
            if line.strip().startswith("c."):
                if line.strip().split(".")[1] not in jupyter_lab_settings:
                    jupyter_lab_settings[line.strip().split(".")[1].strip()] = {}
                try:
                    jupyter_lab_settings[line.strip().split(".")[1]][
                        line.strip().split(".")[2].split("=")[0].strip()
                    ] = eval(line.strip().split("=")[1])
                except:
                    jupyter_lab_settings[line.strip().split(".")[1]][
                        line.strip().split(".")[2].split("=")[0].strip()
                    ] = line.strip().split("=")[1]

    return jupyter_lab_settings
