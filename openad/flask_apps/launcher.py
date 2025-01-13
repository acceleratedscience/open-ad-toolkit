import os
import socket
import webbrowser
import warnings
from time import sleep
from IPython.display import HTML, display
from IPython.display import IFrame as display_iframe
from flask import Flask, send_from_directory
from flask_cors import CORS, cross_origin
from openad.app.global_var_lib import _repo_dir
from openad.helpers.output import output_text, output_error
from openad.helpers.output_msgs import msg
from openad.helpers.general import next_avail_port
from openad.app.global_var_lib import GLOBAL_SETTINGS
import openad.helpers.jupyterlab_settings as jl_settings


def launch(cmd_pointer=None, routes=None, app_name="", query="", hash=""):
    if not routes:
        output_error(msg("err_routes_required"))
        return
    JL_PROXY = False
    FORCE_PROXY = False  # Set this to True to force the use of the proxy for testing (Jupyter only)
    IS_STATIC = ""
    try:
        jl = jl_settings.get_jupyter_lab_config()

        if jl["ServerApp"]["allow_remote_access"] is True and "127.0.0.1" in jl["ServerProxy"]["host_allowlist"]:
            JL_PROXY = True
            IS_STATIC = "/static"
    except Exception as e:
        JL_PROXY = False
        IS_STATIC = ""

    # Initialize Flask app.
    template_folder = os.path.dirname(os.path.abspath(__file__))
    app = Flask("OpenAD", template_folder=template_folder)
    app.config["WTF_CSRF_ENABLED"] = False
    CORS(
        app,
        allow_headers="*",
        origins="*",
        resources={
            r"/api/*": {"origins": "*"},
            r"/js/*": {"origins": "*"},
            r"/assets/*": {"origins": "*"},
            r"/app/*": {"origins": "*"},
        },
    )

    # Make main CSS files available.

    @app.route(f"{IS_STATIC}/css/<path>")
    @cross_origin()
    def static_css(path):
        return send_from_directory(_repo_dir + "/../flask_apps/_css", f"{path}")

    # Make main JS files available.
    @app.route(f"{IS_STATIC}/js/<path>")
    @cross_origin()
    def static_js(path):
        return send_from_directory(_repo_dir + "/../flask_apps/_js", f"{path}")

    @app.route("/css/<path>")
    @cross_origin()
    def css(path):
        return send_from_directory(_repo_dir + "/../flask_apps/_css", f"{path}")

    # Make main JS files available.
    @app.route("/js/<path>")
    @cross_origin()
    def js(path):
        return send_from_directory(_repo_dir + "/../flask_apps/_js", f"{path}")

    # Make app files available.
    flask_dir = os.path.dirname(os.path.abspath(__file__))

    @app.route("/app/<path:subpath>")
    @cross_origin()
    def app_dir(subpath):
        if GLOBAL_SETTINGS["display"] != "notebook":
            suffix = ""
        else:
            suffix = IS_STATIC
        return send_from_directory(f"{flask_dir}/{app_name}{suffix}", f"{subpath}")

    # Unpack routes.
    for route in routes:
        func = routes[route]["func"]
        method = routes[route]["method"] if "method" in routes[route] else "GET"
        cross_origin()(app.route(route, methods=[method])(func))

        # app.route(route, methods=[method])(func)

        # This is the equivalent of:
        # @app.route('/', methods=['GET'])
        # def home():
        #     return render_template('/home.html')

    # Determine port and host.
    host, port = next_avail_port()

    # Launch the UI
    if GLOBAL_SETTINGS["display"] == "notebook":
        # Jupyter --> Render iframe.

        # Rendering the iframe in the traditional way doesn't let us
        # style it, so we have to use a little hack. Jupyter doesn't
        # like our hack, so we have to suppress the warning.
        # The "regular" way of rendering the iframe would be:
        #
        #   from IPython.display import IFrame, display
        #   iframe = IFrame(src=f'http://{host}:{port}', width='100%', height=700)
        #   display(iframe)

        width = "100%"
        height = 700

        prefix = os.environ.get("NB_PREFIX")
        if prefix is None:
            prefix = ""
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=UserWarning)
            if FORCE_PROXY == True or JL_PROXY == True:
                iframe_html = f'<iframe src="{prefix}/proxy/{port}/{query}{hash}" width="{width}" height="{height}" style="border: solid 1px #ddd;"></iframe>'
            else:
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
    output_text(msg("flask_launch", port), pad_top=1)

    # Launch server.
    if GLOBAL_SETTINGS["display"] == "notebook":
        # Jupyter --> Start the Flask app in a separate thread.
        from threading import Thread

        thread = Thread(target=lambda: app.run(host=host, port=port))
        thread.start()
        sleep(0.5)
    else:
        # CLI --> Start the Flask app in the main thread.
        app.run(host=host, port=port)
        return True
