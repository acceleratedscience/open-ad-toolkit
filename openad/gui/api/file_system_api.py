import json
from flask import request
from openad.core.lang_file_system import get_workspace_files as _get_workspace_files
from openad.core.lang_file_system import get_file as _get_file


class FileSystemApi:
    def __init__(self, cmd_pointer):
        self.cmd_pointer = cmd_pointer

    # Fetch the name of the active workspace.
    def get_workspace_name(self):
        return self.cmd_pointer.settings["workspace"].upper()

    # Fetch your active workspace's content as a JSON object.
    def get_workspace_files(self):
        # data = request.data.decode("utf-8")
        data = json.loads(request.data) if request.data else {}
        path = data["path"] if "path" in data else ""
        return _get_workspace_files(self.cmd_pointer, path)

    # Fetch a file's content as a JSON object.
    def get_file(self):
        data = json.loads(request.data) if request.data else {}
        path = data["path"] if "path" in data else ""
        return _get_file(self.cmd_pointer, path)
