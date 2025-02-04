"""
Allow for easy import of Notebook CSS styles.

Note: simply importing the module with `import openad.notebooks.styles`
does not re-initialize the module when re-running the cell. This causes all
the styling to disappear every time you clear all cells output, until you
restart yoru kernel. Importing the module as a Notebook plugin circumvents that,
as load_ipython_extension() is triggered every time the extension is loaded.

Usage:
%reload_ext openad.notebooks.styles
"""

from IPython.display import display
from .init_style import NotebookStyles

nb_styles = NotebookStyles()


def load_ipython_extension(ip):
    display(nb_styles.css())
