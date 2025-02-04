import sys
import importlib
from IPython.core.getipython import get_ipython


def load_ipython_extension(ipython):
    print("SUCCESS")


# Register the extension to be loaded automatically
# ip = get_ipython()
# print(10)
# if ip:
#     print(11)
#     load_ipython_extension(ip)
