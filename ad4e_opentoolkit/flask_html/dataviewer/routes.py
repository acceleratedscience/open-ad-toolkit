from flask import render_template


def fetchRoutes(data_str):

    def home():
        return render_template('/dataviewer/index.html', data=data_str)

    def test():
        return render_template('/dataviewer/test.html')

    routes = {
        '/': {
            'func': home,
            'method': 'GET'
        },
        '/test': {
            'func': test,
            'method': 'GET'
        }
    }
    return routes
