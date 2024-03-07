import os
import json
from flask import request
from openad.workers.file_system import fs_get_workspace_files, fs_get_file
from openad.helpers.output import output_success


class FileSystemApi:
    def __init__(self, cmd_pointer):
        self.cmd_pointer = cmd_pointer

    # Fetch list of workspaces
    def get_workspaces(self):
        return {
            "all": self.cmd_pointer.settings["workspaces"],
            "active": self.cmd_pointer.settings["workspace"],
        }

    # Fetch the name of the active workspace.
    def get_workspace(self):
        return self.cmd_pointer.settings["workspace"]

    # Set the active workspace.
    def set_workspace(self):
        data = json.loads(request.data) if request.data else {}
        new_workspace_name = data["workspace"] if "workspace" in data else ""
        current_workspace_name = self.cmd_pointer.settings["workspace"]
        workspace_diff = new_workspace_name != current_workspace_name
        workspace_exists = new_workspace_name in self.cmd_pointer.settings["workspaces"]
        if workspace_exists and workspace_diff:
            self.cmd_pointer.settings["workspace"] = new_workspace_name
            self.cmd_pointer.histfile = os.path.expanduser(
                self.cmd_pointer.workspace_path(new_workspace_name) + "/.cmd_history"
            )
            output_success(f"Workspace set to <yellow>{new_workspace_name}</yellow>")

        return "ok"

    # Fetch your active workspace's content as a JSON object.
    def get_workspace_files(self):
        # data = request.data.decode("utf-8")
        data = json.loads(request.data) if request.data else {}
        path = data["path"] if "path" in data else ""
        return fs_get_workspace_files(self.cmd_pointer, path)

    # Fetch a file's content as a JSON object.
    def get_file(self):
        data = json.loads(request.data) if request.data else {}
        path = data["path"] if "path" in data else ""
        return fs_get_file(self.cmd_pointer, path)
