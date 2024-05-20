import logging
# import logging.config
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


def get_level(level: int | str):
    """get the correct log level from a string. 'DEBUG'=logging.DEBUG"""
    if isinstance(level.upper(), str):
        return getattr(logging, level.upper())
    else:
        return level


def get_logger(name, level: int = logging.CRITICAL, disable_existing_loggers: bool=False):
    level = get_level(level)
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

    # if disable_existing_loggers:
    #     logging.config.dictConfig({
    #         'version': 1,
    #         'disable_existing_loggers': True,
    #     })
    return logger