from adccl import AD
def load_ipython_extension(ipython):
    ipython.register_magics(AD)
