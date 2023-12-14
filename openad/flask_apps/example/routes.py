import os
import signal
import json
from flask import render_template, request, redirect, url_for
from urllib.parse import parse_qs


def fetchRoutesExample(data):
    # This is the object that is returned when launching the Flask app.
    # When submitting the form, the data is stored here and can be accessed
    # within Jupyter or the CLI.
    result = {}

    #
    #

    def home():
        return render_template("/example/index.html", data=data)

    def foo():
        return render_template("/example/foo.html")

    def bar():
        return render_template("/example/bar.html")

    # Traditional form submit.
    def submit():
        # Fetch the submitte data.
        result["form"] = dict(request.form)

        # Alternative way to get a specific field's data.
        value = request.form.get("someInput")

        # Redirect to success page.
        return redirect(url_for("success", someInput=value))

    # Ajax form submit.
    # - - - - - - - - -
    # Generall this way is preferred because if you use a regular
    # form submit, the url will change causing window.close() to
    # no longer work.
    def submit_ajax():
        # Fetch submitted data.
        result["form"] = json.loads(request.data.decode("utf-8"))
        value = result["form"]["someAjaxInput"]

        # Success message.
        return value

    def success():
        # Kill server
        os.kill(os.getpid(), signal.SIGINT)

        # Fetch submitted value.
        query_string = request.query_string.decode("utf-8")
        query = dict(parse_qs(query_string))
        value = query["someInput"][0] if "someInput" in query else ""

        return render_template("/example/success.html", value=value)

    routes = {
        "/": {
            "func": home,
            "method": "GET",  # Default methos GET doesn't *have* to be specified, just here for clarity.
        },
        "/foo": {"func": foo, "method": "GET"},
        "/bar": {"func": bar, "method": "GET"},
        "/submit": {"func": submit, "method": "POST"},
        "/submit-ajax": {"func": submit_ajax, "method": "POST"},
        "/success": {"func": success, "method": "GET"},
    }

    return routes, result
