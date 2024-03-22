"""
This file contains all the API endpoints consumed by the GUI.

They can be accessed on: http://localhost:8024/api/v1/<endpoint>
"""

from openad.gui.api.general_api import GeneralApi
from openad.gui.api.file_system_api import FileSystemApi
from openad.gui.api.molecules_api import MoleculesApi


def fetchRoutes(cmd_pointer):
    api_v1 = "/api/v1"

    # Import the different API endpoints.
    general_api = GeneralApi(cmd_pointer)
    file_system_api = FileSystemApi(cmd_pointer)
    molecules_api = MoleculesApi(cmd_pointer)

    routes = {
        #
        # General
        f"{api_v1}/": {"func": general_api.landing, "method": "GET"},
        f"{api_v1}/test": {"func": general_api.test, "method": "GET"},
        f"{api_v1}/health": {"func": general_api.health, "method": "GET"},
        f"{api_v1}/exec-command": {"func": general_api.exec_command, "method": "POST"},
        #
        # File system
        f"{api_v1}/get-workspaces": {"func": file_system_api.get_workspaces, "method": "GET"},
        f"{api_v1}/get-workspace": {"func": file_system_api.get_workspace, "method": "GET"},
        f"{api_v1}/set-workspace": {"func": file_system_api.set_workspace, "method": "POST"},
        f"{api_v1}/get-workspace-files": {"func": file_system_api.get_workspace_files, "method": "POST"},
        f"{api_v1}/get-file": {"func": file_system_api.get_file, "method": "POST"},
        #
        # Molecules
        f"{api_v1}/get-mol-data": {"func": molecules_api.get_mol, "method": "POST"},
        f"{api_v1}/get-mol-viz-data": {"func": molecules_api.get_mol_viz_data, "method": "POST"},
        f"{api_v1}/get-molset": {"func": molecules_api.get_molset, "method": "POST"},
        f"{api_v1}/remove-from-molset": {"func": molecules_api.remove_from_molset, "method": "POST"},
        f"{api_v1}/clear-molset-working-copy": {"func": molecules_api.clear_molset_working_copy, "method": "POST"},
        f"{api_v1}/save-molset-changes": {"func": molecules_api.save_molset_changes, "method": "POST"},
        f"{api_v1}/get-mol-data-from-molset": {"func": molecules_api.get_mol_data_from_molset, "method": "POST"},
        #
        # "/": {"func": home, "method": "GET"},
        # "/submit": {"func": submit, "method": "POST"},
        # "/success": {"func": success, "method": "GET"},
    }
    return routes
