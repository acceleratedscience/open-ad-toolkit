# Universal spinner.
# - - -
# Inherits all methods from Halo but sets default parameters and adds some styling.
# Note: text_color='grey' results in black text, so we use our own styling instead.
#
# Usage:
# from openad.helpers.spinner import spinner
# spinner.start("Please hold while we do something...")

from openad.helpers.general import is_notebook_mode
from openad.helpers.output import output_text


if is_notebook_mode():
    from halo import HaloNotebook as Halo
else:
    from halo import Halo


class Spinner(Halo):
    verbose = True

    def __init__(self, verbose=True):
        self.verbose = verbose
        # Alternative spinners:
        # simpleDotsScrolling, interval=100

        if self.verbose is True:
            super().__init__(spinner="triangle", color="white", interval=700)

    def start(self, text=None, no_format=False):
        if self.verbose is True:
            if no_format:
                text = output_text(text, return_val=True, jup_return_format="plain") if text else None
            else:
                text = (
                    output_text(f"<soft>{text}...</soft>", return_val=True, jup_return_format="plain") if text else None
                )
            super().start(text)

    def succeed(self, *args, **kwargs):
        if self.verbose is True:
            return super().succeed(*args, **kwargs)

    def info(self, *args, **kwargs):
        if self.verbose is True:
            super().info(*args, **kwargs)
            return super().start(*args, **kwargs)

    def warn(self, *args, **kwargs):
        if self.verbose is True:
            return super().warn(*args, **kwargs)

    def fail(self, *args, **kwargs):
        if self.verbose is True:
            return super().fail(*args, **kwargs)

    def stop(self):
        if self.verbose is True:
            return super().stop()


spinner = Spinner()
