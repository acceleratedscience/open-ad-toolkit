from ad4e_opentoolkit.flask_html import launcher
from ad4e_opentoolkit.flask_html.molsviewer.routes import fetchRoutes

mol_name_cache = {}  # caches molecule names


def display_mols(cmd_pointer, parser):
    # Load routes and launch browser UI.
    routes, the_mols2grid = fetchRoutes(cmd_pointer, parser)
    if routes and not cmd_pointer.notebook_mode:
        launcher.launch(cmd_pointer, routes, 5000)
    else:
        return the_mols2grid
