from openad.app.magic.openad import AD


def load_ipython_extension(ipython):
    ipython.register_magics(AD)
