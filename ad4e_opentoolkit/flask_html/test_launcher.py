# THIS IS JUST A TEST -- moenen

import threading
import os
import socket
from flask import Flask, render_template, send_from_directory, request
import webbrowser
from multiprocessing import Process

template_folder = os.path.abspath('./')
print(888, template_folder)
app = Flask('Background Flask', template_folder=template_folder)


@app.route('/')
def top():
    return render_template('/test.html')


def run_flask_server():
    app.run(host='127.0.0.1', port=5000, threaded=True)


if __name__ == '__main__':
    socket.setdefaulttimeout(1)

    # Start the Flask server in a separate thread
    flask_thread = threading.Thread(target=run_flask_server)
    flask_thread.daemon = True  # This allows the thread to exit when the main program exits
    flask_thread.start()

    # Open the web browser
    webbrowser.open('http://127.0.0.1:5000', new=1)

    # Keep the main thread running (you can use some form of input here to control when to exit)
    input("Press Enter to exit...")
