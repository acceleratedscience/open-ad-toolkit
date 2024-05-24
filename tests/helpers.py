import random
import string


def random_name(length: int = 10) -> str:
    """generates a random string

    Returns:
        str: random name
    """
    letters = string.ascii_lowercase
    return "pytest_" + "".join(random.choice(letters) for i in range(length))
