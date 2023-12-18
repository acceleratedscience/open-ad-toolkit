from openad.flask_apps import launcher
from openad.flask_apps.example.routes import fetchRoutesExample


def flask_example(cmd_pointer, parser):
    data = {"James": 27, "Amy": 22, "Charlie": 24, "Sarah": 29}

    # Load routes and launch browser UI.
    routes, return_data = fetchRoutesExample(data)

    if cmd_pointer.notebook_mode:
        # Jupyter
        launcher.launch(cmd_pointer, routes, "example")
        return return_data
    else:
        # CLI
        launcher.launch(cmd_pointer, routes, "example")
