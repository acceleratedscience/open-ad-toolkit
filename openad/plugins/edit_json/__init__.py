# Import a class instance instead of the class itself.
import sys
from .edit_json import EditJson

sys.modules[__name__] = EditJson()
