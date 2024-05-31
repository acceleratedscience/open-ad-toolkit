import logging
from collections import OrderedDict
from typing import Generic, Hashable, Optional, TypeVar

# import logging.config
import sys
import os


##########################################################################
# region - Logging
##########################################################################

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
    if isinstance(level, str):
        return getattr(logging, level.upper())
    else:
        return level


def get_logger(
    name,
    level: int = logging.CRITICAL,
    color: bcolors = bcolors.ENDC,
    disable_existing_loggers: bool = False,
):
    level = get_level(level)
    if ENV_DEBUG == "services":
        level = logging.DEBUG
    logger = logging.getLogger(name)
    logger.setLevel(level)
    fmt = (
        "%(asctime)s|"
        + bcolors.OKBLUE
        + "%(levelname)s|%(module)s:%(funcName)s:%(lineno)d| "
        + color
        + "%(message)s"
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


##########################################################################
# region - Lru Cache
##########################################################################

T = TypeVar("T")


class LruCache(Generic[T]):
    def __init__(self, capacity: int):
        self.capacity = capacity
        self.__cache: OrderedDict[Hashable, T] = OrderedDict()

    def get(self, key: Hashable) -> Optional[T]:
        if key not in self.__cache:
            return None
        self.__cache.move_to_end(key)
        return self.__cache[key]

    def insert(self, key: Hashable, value: T) -> None:
        if len(self.__cache) == self.capacity:
            self.__cache.popitem(last=False)
        self.__cache[key] = value
        self.__cache.move_to_end(key)

    def __len__(self) -> int:
        return len(self.__cache)

    def clear(self) -> None:
        self.__cache.clear()
