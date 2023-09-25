import os
import signal
import json
from flask import render_template, request, redirect, url_for


def fetchRoutes(data):

    # This is the object that is returned when launching the Flask app.
    # When submitting the form, the data is stored here and can be accessed
    # within Jupyter or the CLI.
    return_data = {}

    #
    #

    def home():
        return render_template('/example/index.html', data=data)

    def foo():
        return render_template('/example/foo.html')

    def bar():
        return render_template('/example/bar.html')

    def submit():
        # Fetch the submitte data.
        # - - -
        # Note: you can also do this:
        # request.form.get('someInput')
        # return_data = dict(request.form)
        return_data['form'] = dict(request.form)

        # Kill the Flask app.
        os.kill(os.getpid(), signal.SIGINT)

        # Success message.
        return f'You submitted: {return_data["form"]["someInput"]}'

    def success():
        print(1111111)
        return_data['form'] = json.loads(request.data.decode('utf-8'))
        return f'You submitted: {return_data["form"]["someAjaxInput"]}'

    def submit_ajax():
        # return redirect(url_for('success'))

        # Fetch submitted data.
        return_data['form'] = json.loads(request.data.decode('utf-8'))

        # Kill the Flask app.
        os.kill(os.getpid(), signal.SIGINT)

        # Success message.
        return f'You submitted: {return_data["form"]["someAjaxInput"]}'

    routes = {
        '/': {
            'func': home,
            'method': 'GET'  # Default methos GET doesn't need to be specified, just here for clarity.
        },
        '/foo': {
            'func': foo,
            'method': 'GET'
        },
        '/bar': {
            'func': bar,
            'method': 'GET'
        },
        '/submit': {
            'func': submit,
            'method': 'POST'
        },
        '/submit-ajax': {
            'func': submit_ajax,
            'method': 'POST'
        },
        '/success': {
            'func': success,
            'method': 'GET'
        }
    }

    return routes, return_data
