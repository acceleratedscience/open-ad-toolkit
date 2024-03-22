"""
Meant to simplify measuring the duraction of an operation.
Basically because we miss console.time / console.timeEnd from JS.
"""

import time


IDENTIFIERS = {}


def timeit(identifier=None, end=False):
    # End
    if end:
        timeit_end(identifier)

    # Start
    else:
        IDENTIFIERS[identifier] = time.time()


def timeit_end(identifier=None):
    start_time = IDENTIFIERS[identifier]
    del IDENTIFIERS[identifier]

    # Calculate total time
    total_time_ms = round((time.time() - start_time) * 1000)
    if total_time_ms > 1000:
        total_time = f"{total_time_ms / 1000} s"
    else:
        total_time = f"{total_time_ms} ms"

    identifier = f"{identifier} " if identifier else ""
    print(f"Operation {identifier}took {total_time}")
