import os
import socket
import webbrowser
from flask import Flask, render_template, send_from_directory, request
from ad4e_opentoolkit.app.global_var_lib import _repo_dir
from ad4e_opentoolkit.helpers.output import msg, output_text, output_error


def launch(cmd_pointer, routes=None, port=5000, host='127.0.0.1'):
    if not routes:
        output_error('Routes are required to launch Flask server.')
        return

    # Initialize Flask app.
    template_folder = os.path.dirname(os.path.abspath(__file__))
    app = Flask('OpenAD', template_folder=template_folder)

    # Make main CSS files available.
    @app.route('/css/<path>')
    def send_main_css(path):
        return send_from_directory(_repo_dir + '/../flask_apps/_css', f'{path}')

    # Unpack routes.
    for route in routes:
        func = routes[route]['func']
        method = routes[route]['method'] if 'method' in routes[route] else 'GET'
        app.route(route, methods=[method])(func)

        # This is the equivalent of:
        # @app.route('/', methods=['GET'])
        # def home():
        #     return render_template('/home.html')

    # Open browser.
    socket.setdefaulttimeout(1)
    webbrowser.open(f'http://{host}:{port}', new=1)

    # Remove Flask startup message.
    import sys
    cli = sys.modules['flask.cli']
    cli.show_server_banner = lambda *x: None

    # Remove logging of warning & informational messages.
    import logging
    log = logging.getLogger('werkzeug')
    log.setLevel(logging.ERROR)

    # Display our own launch message.
    output_text(msg('flask_launch', 'Data Viewer', port), cmd_pointer, pad_top=1)

    # Launch server.
    app.run(host=host, port=port)
    return True


if __name__ == "__main__":
    launch()
