from openad.app.magic.openad_magic import AD


def load_ipython_extension(ipython):
    ipython.register_magics(AD)
