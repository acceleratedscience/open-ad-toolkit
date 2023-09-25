import os
import socket
import webbrowser
from flask import Flask, send_from_directory
from ad4e_opentoolkit.app.global_var_lib import _repo_dir
from ad4e_opentoolkit.helpers.output import msg, output_text, output_error
from ad4e_opentoolkit.helpers.general import next_avail_port


def launch(cmd_pointer, routes=None):
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
    port, host = next_avail_port()
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
