ABOUT FLASK APPLICATIONS
------------------------

This folder contains all of our flask applications. Use the example app as
a template to hit the ground running when creating new apps. It is launched
from /core/lang_dev --> flask_example()

To launch the example app, run:

    flask example

# Note: Currently, the /submit routes in the example are logging a timeout
# error to the terminal, for which I don't have a good explanation. It doesn't
# affect the rest of the functionality though.
# - moenen Sep 25, 2023

Basic code to launch from python:

    from openad.flask_apps import launcher
    from openad.flask_apps.my_app.routes import fetchRoutesMyApp

    def my_app_command_func(cmd_pointer, parser):
        data = {
            'foo': 123
        }

        # Load routes.
        routes = fetchRoutes(data)

        # Launch browser UI.
        launcher.launch(cmd_pointer, routes, 'my_app')

The launcher script exposes two default routes:
- /css/... -> This exposes the shared _css folder (for carbon + shared styles)
- /app/... -> This exposes the app folder, so you can link to your css and js

Note that in order to expose app, you need to pass the app_name parameter to
the launcher, i.e. the name of the app folder.
