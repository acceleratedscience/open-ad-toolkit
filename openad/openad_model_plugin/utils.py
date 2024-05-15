import logging
import sys
import os


ENV_DEBUG = os.environ.get("debug") or os.environ.get("DEBUG")


class bcolors:
    """Add colors to print statements"""

    HEADER = "\033[95m"
    OKBLUE = "\033[94m"
    OKCYAN = "\033[96m"
    OKGREEN = "\033[92m"
    WARNING = "\033[93m"
    FAIL = "\033[91m"
    ENDC = "\033[0m"
    BOLD = "\033[1m"
    UNDERLINE = "\033[4m"


def get_logger(name, level: int = logging.CRITICAL):
    if ENV_DEBUG:
        level = logging.DEBUG
    logger = logging.getLogger(name)
    logger.setLevel(level)
    fmt = (
        "\n%(asctime)s |"
        + bcolors.OKBLUE
        + " %(levelname)s | %(module)s:%(funcName)s:%(lineno)d |"
        + bcolors.WARNING
        + " %(message)s"
        + bcolors.ENDC
    )
    formatter = logging.Formatter(fmt=fmt, datefmt="%Y/%m/%d %H:%M:%S")
    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    return logger
