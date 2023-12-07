# This file contains snippets of text that are
# repeated throughout the toolkit documentation.
# See toolkit_main.py -> load_toolkit()

NOTE_USING_CLAUSE = """
<bold>Note:</bold> The <cmd>using</cmd> clause requires all enclosed parameters to be defined in the same order as listed below.
"""

SAVE_AS_CLAUSE = """
Use the <cmd>save as</cmd> clause to save the results as a csv file in your current workspace.
"""

HOW_TO_LIST_COLLECTIONS = """
Use the command <cmd>display all collections</cmd> to list available collections.
"""

HOW_TO_LIST_DOMAINS = """
Use the command <cmd>display all collections</cmd> to find available domains.
"""

snippets = {k: v for k, v in globals().items() if not k.startswith("__")}
