"""
This file holds all the GUI-related CLI commands.
"""

from openad.helpers.output import output_text, output_error, output_success, output_warning
from openad.gui.gui_launcher import gui_init


# Command: install gui
def install_gui(cmd_pointer, parser):
    output_error("This is an interaction prototype, installation is not yet supported.")


# Command: launch gui
def launch_gui(cmd_pointer, parser):
    gui_init(cmd_pointer)
