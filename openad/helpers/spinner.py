# Universal spinner.
# - - -
# Inherits all methods from Halo but sets default parameters and adds some styling.
# Note: text_color='grey' results in black text, so we use our own styling instead.

from openad.helpers.general import is_notebook_mode
from openad.helpers.output import msg, output_text

if is_notebook_mode():
    from halo import HaloNotebook as Halo
else:
    from halo import Halo


class Spinner(Halo):
    def __init__(self):
        # Alternative spinners:
        # simpleDotsScrolling, interval=100
        super().__init__(spinner="dots", color="white")

    def start(self, text=None, no_format=False):
        if no_format:
            text = output_text(text, return_val=True, jup_return_format="plain") if text else None
        else:
            text = output_text(f"<soft>{text}...</soft>", return_val=True, jup_return_format="plain") if text else None
        super().start(text)


spinner = Spinner()
