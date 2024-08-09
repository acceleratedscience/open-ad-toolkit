# This file contains snippets of text that are
# repeated throughout the toolkit documentation.
# See toolkit_main.py -> load_toolkit()

REACTION_NOTATION = """
Reactions are defined by combining two SMILES strings delimited by a period. For example: <cmd>'BrBr.c1ccc2cc3ccccc3cc2c1'</cmd>
"""
NOTE_USING_CLAUSE = """
<bold>Note:</bold> The <cmd>using</cmd> clause requires all enclosed parameters to be defined in the same order as listed below.
"""
USING_CLAUSE_OPTIONS = """
Optional Parameters that can be specified in the <cmd>using</cmd> clause:
"""

USING_AI_MODEL = """
<cmd>ai_model='<model_name>'</cmd> What model to use. Use the command <cmd>list rxn models</cmd> to list all available models. The default is '2020-07-01'.
"""

USE_SAVED_CLAUSE = """
You can reuse previously generated results by appending the optional <cmd>use_saved</cmd> clause. This will reuse the results of a previously run command with the same parameters, if available.
"""

snippets = {k: v for k, v in globals().items() if not k.startswith("__")}
