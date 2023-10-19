from flask import render_template


def fetchRoutesDataViewer(data, headers):

    def home():
        return render_template('/dataviewer/index.html', data=data, headers=headers)

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
