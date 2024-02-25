"""
This file contains all the API endpoints consumed by the GUI.

They can be accessed on: http://localhost:5000/api/v1/<endpoint>
"""

from openad.gui.api.general_api import GeneralApi
from openad.gui.api.file_system_api import FileSystemApi
from openad.gui.api.molecules_api import MoleculesApi


def fetchRoutes(cmd_pointer):
    api_v1 = "/api/v1"

    # Import the different API endpoints.
    general = GeneralApi(cmd_pointer)
    file_system = FileSystemApi(cmd_pointer)
    molecules = MoleculesApi(cmd_pointer)

    routes = {
        #
        # General
        f"{api_v1}/test": {"func": general.test, "method": "GET"},
        f"{api_v1}/exec-command": {"func": general.exec_command, "method": "POST"},
        #
        # File system
        f"{api_v1}/get-workspace-name": {"func": file_system.get_workspace_name, "method": "GET"},
        f"{api_v1}/get-workspace-files": {"func": file_system.get_workspace_files, "method": "POST"},
        f"{api_v1}/get-file": {"func": file_system.get_file, "method": "POST"},
        #
        # Molecules
        f"{api_v1}/get-mol-data": {"func": molecules.get_mol_data, "method": "POST"},
        f"{api_v1}/get-mol-viz-data": {"func": molecules.get_mol_viz_data, "method": "POST"},
        #
        # "/": {"func": home, "method": "GET"},
        # "/submit": {"func": submit, "method": "POST"},
        # "/success": {"func": success, "method": "GET"},
    }
    return routes
