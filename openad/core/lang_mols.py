"""automates the fdisplay of mols2grid for the command line and Jupyter"""
from openad.flask_apps import launcher
from openad.flask_apps.molsviewer.routes import fetchRoutesMols2Grid

mol_name_cache = {}  # caches molecule names


def display_mols(cmd_pointer, parser):
    """launches molecule viewer"""
    # Load routes and launch browser UI.
    routes, the_mols2grid = fetchRoutesMols2Grid(cmd_pointer, parser)
    if routes and not cmd_pointer.notebook_mode:
        launcher.launch(cmd_pointer, routes, "molsviewer")
    else:
        return the_mols2grid
