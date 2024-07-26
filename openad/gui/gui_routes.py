"""
This file contains all the API endpoints consumed by the GUI.

They can be accessed on: http://0.0.0.0:8024/api/v1/<endpoint>
"""

from openad.gui.api.general_api import GeneralApi
from openad.gui.api.file_system_api import FileSystemApi
from openad.gui.api.molecules_api import MoleculesApi
from openad.gui.api.result_api import ResultApi
from openad.gui.api.dataframe_api import DataframeApi

# fmt: off
def fetchRoutes(cmd_pointer):
    api_v1 = "/api/v1"

    # Import the different API endpoints.
    general_api = GeneralApi(cmd_pointer)
    file_system_api = FileSystemApi(cmd_pointer)
    molecules_api = MoleculesApi(cmd_pointer)
    result_api = ResultApi(cmd_pointer)
    dataframe_api = DataframeApi(cmd_pointer)


    routes = {
        #
        #
        # General
        f"{api_v1}/": {"func": general_api.landing, "method": "GET"},
        f"{api_v1}/test": {"func": general_api.test, "method": "GET"},
        f"{api_v1}/health": {"func": general_api.health, "method": "GET"},
        f"{api_v1}/exec-command": {"func": general_api.exec_command, "method": "POST"},
        #
        #
        # File system
        f"{api_v1}/get-workspaces": {"func": file_system_api.get_workspaces, "method": "GET"},
        f"{api_v1}/get-workspace": {"func": file_system_api.get_workspace, "method": "GET"},
        f"{api_v1}/set-workspace": {"func": file_system_api.set_workspace, "method": "POST"},
        f"{api_v1}/get-workspace-files": {"func": file_system_api.get_workspace_files, "method": "POST"},
        f"{api_v1}/get-file": {"func": file_system_api.get_file, "method": "POST"},
        f"{api_v1}/open-file-os": {"func": file_system_api.open_file_os, "method": "POST"},
        f"{api_v1}/delete-file": {"func": file_system_api.delete_file, "method": "POST"},
        #
        #
        # Molecules - Mols
        f"{api_v1}/get-mol-data": {"func": molecules_api.get_mol_data, "method": "POST"},
        f"{api_v1}/get-mol-viz-data": {"func": molecules_api.get_mol_viz_data, "method": "POST"},
        f"{api_v1}/get-mol-data-from-molset": {"func": molecules_api.get_mol_data_from_molset, "method": "POST"},
        #
        f"{api_v1}/add-mol-to-mymols": {"func": molecules_api.add_mol_to_list, "method": "POST"},
        f"{api_v1}/remove-mol-from-mymols": {"func": molecules_api.remove_mol_from_list, "method": "POST"},
        f"{api_v1}/check-mol-in-mymols": {"func": molecules_api.check_mol_in_list, "method": "POST"},
        f"{api_v1}/enrich-mol": {"func": molecules_api.enrich_mol, "method": "POST"},
        #
        f"{api_v1}/save-mol-as-json": {"func": molecules_api.save_mol_as_json, "method": "POST"},
        f"{api_v1}/save-mol-as-sdf": {"func": molecules_api.save_mol_as_sdf, "method": "POST"},
        f"{api_v1}/save-mol-as-csv": {"func": molecules_api.save_mol_as_csv, "method": "POST"},
        f"{api_v1}/save-mol-as-mdl": {"func": molecules_api.save_mol_as_mdl, "method": "POST"},
        f"{api_v1}/save-mol-as-smiles": {"func": molecules_api.save_mol_as_smiles, "method": "POST"},
        #
        #
        # Molecules - Molsets
        f"{api_v1}/get-molset": {"func": molecules_api.get_molset, "method": "POST"},
        f"{api_v1}/get-molset-mymols": {"func": molecules_api.get_molset_list, "method": "POST"},
        #
        f"{api_v1}/remove-from-molset": {"func": molecules_api.remove_from_molset, "method": "POST"},
        f"{api_v1}/clear-molset-working-copy": {"func": molecules_api.clear_molset_working_copy, "method": "POST"},
        #
        f"{api_v1}/update-molset": {"func": molecules_api.update_molset, "method": "POST"},
        f"{api_v1}/update-molset-mymols": {"func": molecules_api.update_molset_list, "method": "POST"},
        #
        f"{api_v1}/save-molset-as-json": {"func": molecules_api.save_molset_as_json, "method": "POST"},
        f"{api_v1}/save-molset-as-sdf": {"func": molecules_api.save_molset_as_sdf, "method": "POST"},
        f"{api_v1}/save-molset-as-csv": {"func": molecules_api.save_molset_as_csv, "method": "POST"},
        f"{api_v1}/save-molset-as-smiles": {"func": molecules_api.save_molset_as_smiles, "method": "POST"},
        f"{api_v1}/replace-mol-in-molset": {"func": molecules_api.replace_mol_in_molset, "method": "POST"},
        #
        #
        # Result
        f"{api_v1}/get-result": {"func": result_api.get_result, "method": "POST"},
        f"{api_v1}/update-result-molset": {"func": result_api.update_result_molset, "method": "POST"},
        # f"{api_v1}/update-result-data": {"func": result_api.update_result_data, "method": "POST"}, # placeholder for when dataviewer is integrated.
        #
        #
        # Dataframes
        f"{api_v1}/get-dataframe/<df_name>": {"func": dataframe_api.get_dataframe, "method": "POST"},
        f"{api_v1}/update-dataframe-molset/<df_name>": {"func": dataframe_api.update_dataframe_molset, "method": "POST"},
        # f"{api_v1}/update-dataframe-data/<df_name>": {"func": dataframe_api.update_dataframe_data, "method": "POST"},  # placeholder for when dataviewer is integrated.
        #
        # "/": {"func": home, "method": "GET"},
        # "/submit": {"func": submit, "method": "POST"},
        # "/success": {"func": success, "method": "GET"},
    }
    return routes
