# from .edit_json import EditJson

# edit_json = EditJson()  # create an instance of EditJson

# __all__ = ['edit_json']  # set the default export of the module


import sys
from .edit_json import EditJson

sys.modules[__name__] = EditJson()
