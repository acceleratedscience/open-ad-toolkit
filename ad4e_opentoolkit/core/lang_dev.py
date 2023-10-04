from ad4e_opentoolkit.flask_apps import launcher
from ad4e_opentoolkit.flask_apps.example.routes import fetchRoutes


def flask_example(cmd_pointer, parser):
    data = {
        'James': 27,
        'Amy': 22,
        'Charlie': 24,
        'Sarah': 29
    }

    # Load routes and launch browser UI.
    routes, return_data = fetchRoutes(data)

    if routes and not cmd_pointer.notebook_mode:
        # CLI
        launcher.launch(cmd_pointer, routes)
    else:
        # Jupyter
        launcher.launch(cmd_pointer, routes)
        return return_data
