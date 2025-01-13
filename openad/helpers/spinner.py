"""
Universal spinner
- - -
Inherits all methods from Halo but sets default parameters and adds some styling.
Note: text_color='grey' results in black text, so we use our own styling instead.

Usage:
from openad.helpers.spinner import spinner
spinner.start("Please hold while we do something")
spinner.start("Change the message")
spinner.succeed("Done")
spinner.fail("Done")
spinner.stop()
"""

from time import sleep
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

        # Fancy spinner, but requires more CPU, blocking the main thread
        # To do: see is separating thread for spinner resolves this
        wave_spinner = {
            "interval": 700,
            "frames": [
                "▉▋▍▎▏▏",
                "▉▉▋▍▎▏",
                "▋▉▉▋▍▎",
                "▍▋▉▉▋▍",
                "▏▎▋▉▉▋",
                "▏▎▍▋▉▉",
                "▎▏▎▍▋▉",
                "▍▎▏▎▍▋",
                "▋▍▎▏▎▍",
            ],
        }

        if self.verbose is True:
            super().__init__(spinner="triangle", color="white", interval=700)

            # Fancy spinner
            # super().__init__(spinner=wave_spinner, color="yellow")

            # Alternative spinners:
            # https://github.com/sindresorhus/cli-spinners/blob/dac4fc6571059bb9e9bc204711e9dfe8f72e5c6f/spinners.json

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

    def countdown(
        self,
        seconds: int,
        msg: str = None,
        stop_msg: str = None,
    ) -> bool:
        """
        Spinner with countdown timer.

        Parameters
        ----------
        seconds : int
            Number of seconds to countdown from.
        msg : str, optional
            Message to display, with {sec} as placeholder for seconds.
        stop_msg : str, optional
            Message to display when countdown is complete,
            instead of stopping spinner.
        """

        msg = msg or "Waiting {sec} seconds before retrying"
        spinner.start(msg.format(sec=seconds))
        sleep(1)
        if seconds > 1:
            self.countdown(seconds - 1, msg, stop_msg)
        else:
            if stop_msg:
                spinner.start(stop_msg)
            else:
                spinner.stop()
            return True


spinner = Spinner()
