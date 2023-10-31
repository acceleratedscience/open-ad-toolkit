from ad4e_opentoolkit.flask_apps import launcher
from ad4e_opentoolkit.flask_apps.molsgrid.routes import fetchRoutesMols2Grid

mol_name_cache = {}  # caches molecule names


def display_mols(cmd_pointer, parser):
    # Load routes and launch browser UI.
    routes, the_mols2grid = fetchRoutesMols2Grid(cmd_pointer, parser)
    if routes and not cmd_pointer.notebook_mode:
        launcher.launch(cmd_pointer, routes, "molsgrid")
    else:
        return the_mols2grid
