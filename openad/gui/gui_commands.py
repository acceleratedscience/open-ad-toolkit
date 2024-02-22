"""
This file holds all the GUI-related CLI commands.
"""

from openad.helpers.output import output_text, output_error, output_success, output_warning
from openad.gui.gui_launcher import gui_init, gui_shutdown


# Install gui
def install_gui(cmd_pointer, parser):
    output_error("This is an interaction prototype, installation is not yet supported.")


# Launch general gui
def launch_gui(cmd_pointer, parser):
    gui_init(cmd_pointer)


# Launch specific module
def launch_gui_module(cmd_pointer, parser):
    module_name = parser["module_name"]
    if module_name == "filebrowser":
        gui_init(cmd_pointer, "filebrowser")  # Refers to the vue template's filename
    elif module_name == "molviewer":
        gui_init(cmd_pointer, "molviewer")  # Refers to the vue template's filename


# restart the gui server.
def restart_gui(cmd_pointer, parser):
    gui_shutdown()
    gui_init(cmd_pointer, silent=True)


# Terminate the gui server.
def quit_gui(cmd_pointer, parser):
    gui_shutdown()
