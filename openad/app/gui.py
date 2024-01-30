import os
import sys
import html
import socket
import logging
import webbrowser
from pathlib import Path
from threading import Thread
from werkzeug.serving import make_server
from flask import Flask, send_from_directory
from openad.helpers.output_msgs import msg
from openad.helpers.output import output_text, output_error, output_success
from openad.helpers.general import next_avail_port, confirm_prompt
from openad.app.global_var_lib import GLOBAL_SETTINGS

GUI_SERVER = None


def init_gui(cmd_pointer=None):
    """
    Check if the GUI is installed and start the server.

    If this GUI is not installed (openad-gui folder not present),
    suggest to install it, unless the user has chosen before not
    to be asked again (openad-gui folder present without index.html).
    """
    template_folder = Path(__file__).resolve().parents[2] / "openad-gui"

    # GUI is installed, launch it.
    if os.path.exists(template_folder):
        is_installed = template_folder / "index.html"
        if is_installed:
            _launch(cmd_pointer)

    # GUI is not yet installed, suggest installation.
    else:
        install_now = confirm_prompt(
            "The OpenAD GUI (graphical user interface) is not yet installed. Would you like to install it now?"
        )
        if install_now:
            output_error("This is an interaction prototype, installation is not yet supported.")
        else:
            output_text("You can install the GUI at any time by running <cmd>install gui</cmd>")
            remind_me = confirm_prompt("Remind you next time?", default=True)
            if not remind_me:
                # Create the openad-cli folder with an empty README.txt file inside.
                try:
                    Path(template_folder).mkdir(parents=False)
                    readme = Path(template_folder) / "README.txt"
                    readme.touch()
                    readme_text = "You have declined to install the OpenAD GUI.\nYou can install it at any time by running `install gui`.\n- - -\nThe graphical user interface allows you to browse your workspace and visualize your datasets and molecules."
                    readme.write_text(readme_text)
                except BaseException as err:
                    output_error(["Something went wrong while creating the openad-gui folder.", err])


def _launch(cmd_pointer=None, query="", hash=""):
    """
    Launch the GUI web server in a separate thread.
    """

    # if not routes:
    #     output_error(msg("err_routes_required"))
    #     return

    # Initialize Flask app.
    template_folder = Path(__file__).resolve().parents[2] / "openad-gui"
    if not template_folder.exists():
        output_error("The OpenAD GUI folder is missing")
        return
    app = Flask("OpenAD", template_folder=template_folder)

    # Make asset files available (CSS/JS).
    @app.route("/assets/<path>")
    def assets(path):
        return send_from_directory(template_folder / "assets", f"{path}")

    # Shutdown path.
    @app.route("/shutdown", methods=["GET"])
    def shut_down():
        # Shutdown here needs to happen with a delay,
        # because this API call needs to finish before
        # the shutdown will execute.
        def delayed_shutdown():
            import time

            time.sleep(1)
            global GUI_SERVER
            GUI_SERVER.shutdown()

        Thread(target=delayed_shutdown).start()
        return "OpenAD GUI shutdown complete"

    # Serve all other paths.
    @app.route("/", methods=["GET"])
    @app.route("/<path:path>")
    def serve(path=""):
        if path != "" and (template_folder / path).exists():
            return send_from_directory(template_folder, path)
        return send_from_directory(template_folder, "index.html")

    # Determine port and host.
    host, port = next_avail_port()

    # Launch the UI
    if GLOBAL_SETTINGS["display"] == "notebook":
        # Jupyter --> Render iframe.

        # Rendering the iframe in the traditional way doesn't let
        # us style it, so we have to use a little hack, rendering
        # our iframe using HTML. Jupyter doesn't like our hack, so
        # we also have to suppress the warning.
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
    cli = sys.modules["flask.cli"]
    cli.show_server_banner = lambda *x: None

    # Remove logging of warning & informational messages.
    log = logging.getLogger("werkzeug")
    log.setLevel(logging.ERROR)

    # Display our own launch message.
    launch_msg = _render_launch_msg(host, port)
    output_text(launch_msg, pad=1, tabs=1)

    # We need to spin up the persistent GUI server in a separate
    # thread, otherwise it would be blocking the application.
    #
    # Hence, we can't do:
    #   app.run(host=host, port=port)
    #
    # To spin it up in a separate thread, we could do:
    #   THREAD = Thread(target=lambda: app.run(host=host, port=port))
    #
    # But unfortunately Flask's built-in server doesn't provide an option
    # to shut it down programatically from another thread, so instead
    # we use Werkzeug's more advanced WSGI server.
    global GUI_SERVER
    GUI_SERVER = ServerThread(app, host, port)
    GUI_SERVER.start()


# This class replaces Flask's builtin app.run() server,
# and allows us to shut down the server from another thread.
# - - -
# In practice, this happens on KeyboardInterrupt in cmd_line().
# See GUI_SERVER.shutdown() in main.py.
class ServerThread(Thread):
    def __init__(self, app, host, port):
        Thread.__init__(self)
        self.srv = make_server(host, port, app)
        self.ctx = app.app_context()
        self.ctx.push()
        self.host = host
        self.port = port

    def run(self):
        # print("Server started")
        self.srv.serve_forever()

    def shutdown(self):
        self.srv.shutdown()  # Shutdown server
        self.join()  # Close thread
        prefix = f"<red>{html.unescape('&empty;')}</red> "
        output_success([f"{prefix}OpenAD GUI shutdown complete", f"{prefix}{self.host}:{self.port}"], tabs=1)


# Generated a stylized launch message for the web server.
def _render_launch_msg(host, port):
    text = "User interface:"
    url = f"{host}:{port}"
    width = max(len(text), len(url))
    edge = width * "-"
    launch_msg = [
        f"+ {edge} +",
        f"| <yellow>{text}{(width - len(text)) * ' '}</yellow> |",
        f"| <link>{url}</link>{(width - len(url)) * ' '} |",
        f"+ {edge} +",
    ]
    launch_msg = "\n".join(launch_msg)
    return launch_msg


if __name__ == "__main__":
    init_gui()
