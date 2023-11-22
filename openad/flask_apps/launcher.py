import os
import socket
import webbrowser
from flask import Flask, send_from_directory
from openad.app.global_var_lib import _repo_dir
from openad.helpers.output import msg, output_text, output_error
from openad.helpers.general import next_avail_port


def launch(cmd_pointer=None, routes=None, app_name="", query="", hash=""):
    if not routes:
        output_error("Routes are required to launch Flask server.")
        return

    # Initialize Flask app.
    template_folder = os.path.dirname(os.path.abspath(__file__))
    app = Flask("OpenAD", template_folder=template_folder)

    # Make main CSS files available.
    @app.route("/css/<path>")
    def css(path):
        return send_from_directory(_repo_dir + "/../flask_apps/_css", f"{path}")

    # Make main JS files available.
    @app.route("/js/<path>")
    def js(path):
        return send_from_directory(_repo_dir + "/../flask_apps/_js", f"{path}")

    # Make app files available.
    flask_dir = os.path.dirname(os.path.abspath(__file__))

    @app.route("/app/<path:subpath>")
    def app_dir(subpath):
        return send_from_directory(f"{flask_dir}/{app_name}", f"{subpath}")

    # Unpack routes.
    for route in routes:
        func = routes[route]["func"]
        method = routes[route]["method"] if "method" in routes[route] else "GET"
        app.route(route, methods=[method])(func)

        # This is the equivalent of:
        # @app.route('/', methods=['GET'])
        # def home():
        #     return render_template('/home.html')

    # Determine port and host.
    port, host = next_avail_port()

    # Launch the UI
    if cmd_pointer.notebook_mode:
        # Jupyter --> Render iframe.

        # Rendering the iframe in the traditional way doesn't let us
        # style it, so we have to use a little hack. Jupyter doesn't
        # like our hack, so we have to suppress the warning.
        # The "regular" way of rendering the iframe would be:
        #
        #   from IPython.display import IFrame, display
        #   iframe = IFrame(src=f'http://{host}:{port}', width='100%', height=700)
        #   display(iframe)

        import warnings
        from IPython.display import HTML, display

        width = "100%"
        height = 700

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=UserWarning)
            iframe_html = f'<iframe src="http://{host}:{port}{query}{hash}" width="{width}" height="{height}" style="border: solid 1px #ddd;"></iframe>'
            display(HTML(iframe_html))
    else:
        # CLI --> Open browser.
        socket.setdefaulttimeout(1)
        webbrowser.open(f"http://{host}:{port}{query}{hash}", new=1)

    # Remove Flask startup message.
    import sys

    cli = sys.modules["flask.cli"]
    cli.show_server_banner = lambda *x: None

    # Remove logging of warning & informational messages.
    import logging

    log = logging.getLogger("werkzeug")
    log.setLevel(logging.ERROR)

    # Display our own launch message.
    output_text(msg("flask_launch", "Data Viewer", port), cmd_pointer, pad_top=1)

    # Launch server.
    if cmd_pointer.notebook_mode:
        # Jupyter --> Start the Flask app in a separate thread.
        from threading import Thread

        thread = Thread(target=lambda: app.run(host=host, port=port))
        thread.start()
    else:
        # CLI --> Start the Flask app in the main thread.
        app.run(host=host, port=port)
        return True
