import os
import socket
from flask import Flask, render_template, send_from_directory, request
from ad4e_opentoolkit.app.global_var_lib import _repo_dir
from ad4e_opentoolkit.helpers.output import output_error
import webbrowser
from multiprocessing import Process


def launch(routes=None):
    if not routes:
        output_error('Routes are required to launch Flask server.')
        return

    # Initialize Flask app.
    template_folder = os.path.dirname(os.path.abspath(__file__))
    print(444, template_folder)
    app = Flask('OpenAD', template_folder=template_folder)

    # Make main CSS files available.
    @app.route('/css/<path>')
    def send_main_css(path):
        return send_from_directory(_repo_dir + '/../flask_html/css', f'{path}')

    # Unpack routes.
    for route in routes:
        func = routes[route]
        app.route(route)(func)

        # @app.route(route)
        # def dynamic_route():
        #     return func()

    # @app.route('/')
    # def top():
    #     return render_template('/test.html')

    # Open browser.
    socket.setdefaulttimeout(1)
    webbrowser.open('http://127.0.0.1:5000', new=1)

    # Remove logging of warning & informational messages.
    import logging
    log = logging.getLogger('werkzeug')
    log.setLevel(logging.ERROR)

    # Launch server.
    app.run()
    return True


if __name__ == "__main__":
    launch()
