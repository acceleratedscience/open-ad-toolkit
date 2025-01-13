"""
Localization helper functions. This allows us to display text in a user's preferred language.


See your current locale settings:
    run `locale` in your terminal

See available locales:
    run `locale -a` in your terminal

Change session locale to test:
    run `LC_ALL=fr_FR.UTF-8` in your terminal
"""

import os
import locale
from openad.helpers.files import open_file

# https://docs.python.org/3/library/locale.html#locale.setlocale
locale.setlocale(locale.LC_ALL, "")


def localize(text: str | dict):
    """
    Return the relevant localization from a dictionary with localized strings.

    Parameters
    ----------
    Input: str | dict
        foo: {
            "en": "This is English",
            "fr": "Ceci est français général",
            "fr_BE": "Ceci est français belge",
        }
    """
    if not text:
        return None

    # Simple string text
    elif isinstance(text, str):
        return text

    # Localized text
    elif isinstance(text, dict):
        lang = _get_locale("lang")
        region = _get_locale("region")
        _note = None
        if lang:
            _note = text.get(f"{lang}_{region}", None) if region else None
            if not _note:
                _note = text.get(lang, None)
        if not _note:
            _note = text.get("en", None)

        return _note

    return text


def localize_path(path: str, include_region=False):
    """
    Return the path to a localized txt file.

    Parameters
    ----------
    path: str
        Path to the txt file, eg. "path/to/file.txt"
    include_region: bool, optional
        Include the region in the returned path if available.
        This will return "file_en_GB.txt" instead of "file_en.txt"

    Returns
    -------
    The localized path to the txt file, eg. "path/to/file_en_GB.txt"

    """
    lang = _get_locale("lang")
    if lang:
        region = _get_locale("region") if include_region else None
        if region:
            path = path.replace(".txt", f"_{lang}_{region}.txt")
        else:
            path = path.replace(".txt", f"_{lang}.txt")
    return path


def open_localized_file(path: str):
    """
    Look for localized version of a txt file based on
    the user's locale settings and return its content.

    Parameters
    ----------
    path: str
        Path to the txt file, eg. "path/to/file.txt"

    Returns
    -------
    By priority, is it exists:
    1. Language + region, eg. "file_en_GB.txt"
    2. Language, eg. "file_en.txt"
    3. Default, eg. "file.txt"
    """
    if not path:
        return None
    path_localized = localize_path(path)
    path_localized_region = localize_path(path, True)

    # Check if language+region-localized description file exist
    if os.path.isfile(path_localized_region):
        path = path_localized_region

    # Check if language-localized description file exist
    elif os.path.isfile(path_localized):
        path = path_localized

    description = open_file(path)
    return description


def _get_locale(key=None):
    """
    Find out the locale setting of your terminal.

    Parameters
    ----------
    key: str, optional
        Specify which part of the locale to return: "lang", "region", "encoding"
        If no key is specified, the full locale info is returned as a tuple.

    Returns:
    - None if no locale is set.
    - (<language>, <region>, <encoding>) if no key is specified
    - <language> if key="lang"
    - <region> if key="region"
    - <encoding> if key="encoding"

    Locale breakdown:
    fr_FR.ISO8859-15 --> ('fr', 'FR', 'ISO8859-15')
    fr_BE.ISO8859-1 --> ('fr', 'BE', 'ISO8859-1')
    fr_BE.UTF-8 --> ('fr', 'BE', 'UTF-8')
    fr_BE --> ('fr', 'BE', 'ISO8859-1') # Returns default encoding when not specified
    """

    output = None

    try:
        lang = locale.getdefaultlocale()
        output = {
            "lang": lang[0].split("_")[0],
            "region": lang[0].split("_")[1],
            "encoding": lang[1],
        }
    except Exception:  # pylint: disable=broad-except
        pass

    if output:
        if key:
            return output[key]
        else:
            return (output["lang"], output["region"], output["encoding"])
    else:
        return None
