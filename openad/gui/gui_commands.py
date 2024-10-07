"""
This file holds all the GUI-related CLI commands.
"""

from openad.helpers.output import output_text, output_error, output_success, output_warning
from openad.gui.gui_launcher import gui_init, gui_shutdown


# Install gui
def install_gui(cmd_pointer, parser):
    output_error("This is an interaction prototype, installation is not yet supported.")


# Launch module
def launch_gui(cmd_pointer, parser):
    path = parser["path"] if "path" in parser else "~/"
    if path == "filebrowser":
        path = "~/"
    gui_init(cmd_pointer, path)


# restart the gui server.
def restart_gui(cmd_pointer, parser):
    gui_shutdown()
    gui_init(cmd_pointer, silent=True)


# Terminate the gui server.
def quit_gui(cmd_pointer, parser):
    # gui_shutdown()
    import openad.gui.ws_server as ws_server
    import asyncio

    asyncio.run(ws_server.send("HOLY COW"))
    print(">>")
    # print(">> holy COW!")
