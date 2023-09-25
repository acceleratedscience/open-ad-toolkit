from flask import render_template, request


def fetchRoutes(data_str):

    def home():
        return render_template('/example/index.html', data=data_str)

    def foo():
        return render_template('/example/foo.html')

    def bar():
        return render_template('/example/bar.html')

    def submit():
        # Fetch the submitte data.
        data = dict(request.form)

        # Do some stuff here.
        return f'You submitted: {data["someInput"]}'

    def submit_ajax():
        return 'success'

    routes = {
        '/': {
            'func': home,
            'method': 'GET'
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
        }
    }
    return routes
