"""Builds the grammar for the DSL"""

# Note: this file is organized in regions, which collapse in VS Code:
# - Collapse all regions:       hold cmd, then hit K followed by 1
# - Expand all regions:         hold cmd, then hit K followed by J
# - Collapse to any level:      first expand everything, then hold cmd, then hit K followed by any number

import os, re
import glob


# Globals
from pyparsing import (
    Word,
    delimitedList,
    alphas,
    alphanums,
    OneOrMore,
    ZeroOrMore,
    CharsNotIn,
    Forward,
    CaselessKeyword,
    QuotedString,
    ParserElement,
    Suppress,
    Optional,
    Group,
    nums,
    printables,
    # Literal,
    # replaceWith,
    # Combine,
    # pyparsing_test,
    # ParseException,
)

# Main
from openad.core.help import help_dict_create, organize_commands
import openad.toolkit.toolkit_main as toolkit_main  # Not using "from" to avoid circular import.
from openad.smols.smol_grammar import smol_grammar_add
from openad.mmols.mmol_grammar import mmol_grammar_add
from openad.plugins.style_parser import tags_to_markdown


# Helpers
from openad.helpers.general import is_notebook_mode
from openad.helpers.output import output_error, output_text
from openad.helpers.output_msgs import msg
from openad.openad_model_plugin.openad_model_toolkit import service_grammar_add

from openad.openad_model_plugin.catalog_model_services import get_cataloged_service_defs, service_catalog_grammar


# Global variables
from openad.app.global_var_lib import _all_toolkits


(
    get,
    lister,
    description,
    using,
    create,
    s_et,
    unset,
    workspace,
    workspaces,
    context,
    jobs,
    e_xec,
    a_s,
    optimize,
    w_ith,
    toolkits,
    toolkit,
    GPU,
    experiment,
    add,
    run,
    save,
    runs,
    show,
    o_pen,
    mol,
    molecules,
    file,
    d_isplay,
    history,
    data,
    remove,
    update,
    result,
    install,
    launch,
    restart,
    q_uit,
    gui,
    filebrowser,
    molviewer,
) = map(
    CaselessKeyword,
    "get list description using create set unset workspace workspaces context jobs exec\
    as optimize with toolkits toolkit gpu experiment add run save runs show open mol molecules\
    file display history data remove update result install launch restart quit gui filebrowser molviewer".split(),
)
STRING_VALUE = alphanums

##########################################################################

skobki = "(" + ZeroOrMore(CharsNotIn(")")) + ")"
field_def = OneOrMore(Word(alphas, alphanums + "_\"':-") | skobki)


def field_act(s, loc, tok):
    """formatting for fields"""
    return ("".join(tok)).replace('"', '\\"')
    # return ("<" + tok[0] + "> " + " ".join(tok)).replace('"', '\\"')


field_def.setParseAction(field_act)
field_list_def = delimitedList(field_def)


def field_list_act(toks):
    """join toklens with ,"""
    return " , ".join(toks)


field_list_def.setParseAction(field_list_act)

desc = QuotedString("'", escQuote="\\")


ParserElement.inlineLiteralsUsing(Suppress)


name_expr = Word(alphanums + "_" + ".")
key_val_expr = Word(alphanums + "_" + ".")
key_val_expr_num = Word(nums)
key_val_expr_alpha = Word(alphanums + "_" + ".")

key_val_line = Group(name_expr("key") + Suppress("=") + key_val_expr("val"))
key_val_lines = OneOrMore(key_val_line)
complex_obj = Forward()
complex_obj << Group(name_expr("obj_name") + "{" + ZeroOrMore(key_val_lines | complex_obj) + "}")("objects")
GPU_CLAUSE = Suppress(optimize + w_ith) + GPU
AS_EXPERIMENT = Suppress(a_s + experiment) + Word(alphas, alphanums + "_" + ".")
obj = Forward()
obj << Suppress("(") + ZeroOrMore(key_val_lines) + Suppress(")")


delim_value = Group(delimitedList(STRING_VALUE))("list")

##########################################################################

statement = Forward()

statements = []  # Statement definitions
grammar_help = []  # Help text


# Add small molecule grammar
# Note: moved to create_statements()
# smol_grammar_add(statements=statements, grammar_help=grammar_help)

# Add macromolecule grammar
mmol_grammar_add(statements=statements, grammar_help=grammar_help)

# TODO: Organize all other grammar also by individual files

##########################################################################
# region - General
##########################################################################

# General splash screen
statements.append(Forward(CaselessKeyword("openad"))("welcome"))
grammar_help.append(
    help_dict_create(
        name="welcome", category="General", command="openad", description="Display the openad splash screen."
    )
)

# Get status
statements.append(Forward(get + CaselessKeyword("status"))("get_status"))
grammar_help.append(
    help_dict_create(
        name="get status",
        category="General",
        command="get status",
        description="Display the currently selected workspace and toolkit.",
    )
)

# Display history
statements.append(Forward(d_isplay + history("history"))("display_history"))
grammar_help.append(
    help_dict_create(
        name="Display History",
        category="General",
        command="display history ",
        description="Display the last 30 commands run in your current workspace.",
    )
)

# Clear sessions
statements.append(Forward(CaselessKeyword("clear") + CaselessKeyword("sessions"))("clear_sessions"))
grammar_help.append(
    help_dict_create(
        name="Clear Sessions",
        category="General",
        command="clear sessions",
        description="Clear any other sessions that may be running.",
    )
)

# endregion

##########################################################################
# region - Workspaces
##########################################################################

INFO_WORKSPACES = "<soft>To learn more about workspaces, run <cmd>workspace ?</cmd></soft>"

# Set workspaces
statements.append(
    Forward(s_et + workspace("workspace") + Word(alphas, alphanums + "_")("Workspace_Name"))("set_workspace_statement")
)
grammar_help.append(
    help_dict_create(
        name="set workspace",
        category="Workspaces",
        command="set workspace <workspace_name>",
        description=f"Change the current workspace.",
        note=INFO_WORKSPACES,
    )
)

# Get workspace
statements.append(
    Forward(get + workspace("workspace") + Optional(Word(alphas, alphanums + "_")("Workspace_Name")))("get_workspace")
)
grammar_help.append(
    help_dict_create(
        name="get workspace",
        category="Workspaces",
        command="get workspace [ <workspace_name> ]",
        description="Display details a workspace. When no workspace name is passed, details of your current workspace are displayed.",
        note=INFO_WORKSPACES,
    )
)

# Create workspaces
statements.append(
    Forward(
        create
        + workspace("workspace")
        + Word(alphas, alphanums + "_")("Workspace_Name")
        + Optional(description("description") + Suppress("(") + desc("proj_desc") + Suppress(")"))("description_option")
        + Optional(CaselessKeyword("on") + CaselessKeyword("path") + desc("w_path"))("workspace_path")
    )("create_workspace_statement")
)
grammar_help.append(
    help_dict_create(
        name="create workspace",
        category="Workspaces",
        command="create workspace <workspace_name> [ description('<description>') on path '<path>' ]",
        description="Create a new workspace with an optional description and path.",
        note=INFO_WORKSPACES,
    )
)

# Remove workspaces
statements.append(
    Forward(remove + workspace("workspace") + Word(alphas, alphanums + "_")("Workspace_Name"))(
        "remove_workspace_statement"
    )
)
grammar_help.append(
    help_dict_create(
        name="remove workspace",
        category="Workspaces",
        command="remove workspace <workspace_name> ",
        description="Remove a workspace from your registry. Note that this doesn't remove the workspace's directory.",
        note=INFO_WORKSPACES,
    )
)

# list workspaces
statements.append(Forward(lister + workspaces("workspaces"))("list_workspaces"))
grammar_help.append(
    help_dict_create(
        name="list workspaces",
        category="Workspaces",
        command="list workspaces",
        description="Lists all your workspaces.",
        note=INFO_WORKSPACES,
    )
)

# endregion

##########################################################################
# region - Toolkits
# Note toolkits is the Caseless key word now .. simply changed in metadata
##########################################################################

NOTE_TOOLKITS_SEE_ALL = "<soft>To see all available toolkits, run <cmd>list all toolkits</cmd>.</soft>"
NOTE_TOOLKITS = "<soft>To learn more about toolkits, run <cmd>toolkit ?</cmd>.</soft>"


# Toolkit overview screens (explains what, how to install, etc.)
for tk in _all_toolkits:
    # Note: statement creation is moved to create_statements()
    # because they need to be appended after the toolkit command
    # statements, otherwise `<prefix> abc` would be caught by
    # the toolkit command `<prefix>` and cause error.
    # statements.append(Forward(CaselessKeyword(tk))(tk))
    grammar_help.append(
        help_dict_create(
            name=f"{tk} splash",
            category="Toolkits",
            command=tk.lower(),
            # This description is never read. Inside main.py -> do_help()
            # there is a clause that intercepts this command and displays
            # the available commands for the toolkit instead.
            description=f"Display the splash screen for the {tk} toolkit.",
        )
    )

# List toolkits
statements.append(Forward(lister + toolkits("toolkits"))("list_toolkits"))
grammar_help.append(
    help_dict_create(
        name="list toolkits",
        category="Toolkits",
        command="list toolkits",
        description=f"List all installed toolkits.",
        note=f"{NOTE_TOOLKITS_SEE_ALL}\n{NOTE_TOOLKITS}",
    )
)

# List all toolkits
statements.append(Forward(lister + CaselessKeyword("all") + toolkits("toolkits"))("list_all_toolkits"))
grammar_help.append(
    help_dict_create(
        name="list all toolkits",
        category="Toolkits",
        command="list all toolkits",
        description="List all available toolkits.",
        note=NOTE_TOOLKITS,
    )
)

# Install toolkit
# Note: Update toolkit yet to be implemented (currently de-register, then add toolkit with new source)
statements.append(Forward(add + toolkit("workspace") + Word(alphas, alphanums + "_")("toolkit_name"))("add_toolkit"))
grammar_help.append(
    help_dict_create(
        name="add toolkit",
        category="Toolkits",
        command="add toolkit <toolkit_name>",
        description="Install a toolkit.",
        note=f"{NOTE_TOOLKITS_SEE_ALL}\n{NOTE_TOOLKITS}",
    )
)

# Remove toolkit
statements.append(Forward(remove + toolkit + Word(alphas, alphanums + "_")("toolkit_name"))("remove_toolkit"))
grammar_help.append(
    help_dict_create(
        name="remove toolkit",
        category="Toolkits",
        command="remove toolkit <toolkit_name>",
        description=(
            "Remove a toolkit from the registry.\n\n"
            "<b>Note:</b> This doesn't delete the toolkit code. If the toolkit is added again, a backup of the previous install is created in the toolkit directory at <cmd>~/.openad/toolkits</cmd>."
        ),
        note=NOTE_TOOLKITS,
        # Correct description but we have to update the functionality first.
        # description="Remove a toolkit from the registry. This affects all workspaces. A backup of the toolkit directory is stored in <yellow>~/.openad/toolkits_archive</yellow>."
    )
)

# Update toolkit
statements.append(Forward(update + toolkit + Word(alphas, alphanums + "_")("toolkit_name"))("update_toolkit"))
grammar_help.append(
    help_dict_create(
        name="update toolkit",
        category="Toolkits",
        command="update toolkit <toolkit_name>",
        description=("Update a toolkit with the latest version. It is recommended to do this on a regular basis."),
        note=NOTE_TOOLKITS,
    )
)

# Update all toolkits
statements.append(Forward(update + CaselessKeyword("all") + toolkits("toolkits"))("update_all_toolkits"))
grammar_help.append(
    help_dict_create(
        name="update all toolkits",
        category="Toolkits",
        command="update all toolkits",
        description=(
            "Update all installed toolkits with the latest version. Happens automatically whenever OpenAD is updated to a new version."
        ),
        note=NOTE_TOOLKITS,
    )
)

# Set a toolkit as the current context
statements.append(
    Forward(
        s_et
        + context("context")
        + Word(alphas, alphanums + "_")("toolkit_name")
        + Optional(CaselessKeyword("reset"))("reset")
    )("set_context")
)
grammar_help.append(
    help_dict_create(
        name="set context",
        category="Toolkits",
        command="set context <toolkit_name> [ reset ]",
        description="Set your context to the chosen toolkit. By setting the context, the selected toolkit functions become available to you. The optional parameter <cmd>reset</cmd> can be used to reset your login information.",
        note=NOTE_TOOLKITS,
    )
)

# Get the current context
statements.append(Forward(get + context("context"))("get_context"))
grammar_help.append(
    help_dict_create(
        name="get context",
        category="Toolkits",
        command="get context",
        description="Display the currently selected toolkit.",
    )
)

# Unset a toolkit as the current context
statements.append(Forward(unset + context("context"))("unset_context"))
grammar_help.append(
    help_dict_create(
        name="unset context",
        category="Toolkits",
        command="unset context",
        description="Exit your toolkit context. You will no longer have access to toolkit-specific functions.",
        note=NOTE_TOOLKITS,
    )
)

# endregion

##########################################################################
# region - Runs
##########################################################################

NOTE_RUNS = "<soft>To learn more about runs, run <cmd>run ?</cmd>.</soft>"

# Create run
statements.append(Forward(create + run("run"))("create_run"))
grammar_help.append(
    help_dict_create(
        name="create run",
        category="Runs",
        command="create run",
        description="Start recording a run.",
        note=NOTE_RUNS,
    )
)


# Remove run
statements.append(Forward(remove + run("run") + Word(alphas, alphanums + "_")("run_name"))("remove_run"))
grammar_help.append(
    help_dict_create(
        name="remove run",
        category="Runs",
        command="remove run <run_name>",
        description="remove a run.",
        note=NOTE_RUNS,
    )
)


# Save run
statements.append(Forward(save + run("run") + a_s + Word(alphas, alphanums + "_")("run_name"))("save_run"))
grammar_help.append(
    help_dict_create(
        name="save run",
        category="Runs",
        command="save run as <run_name>",
        description="Stop recording a run and save it.",
        note=NOTE_RUNS,
    )
)

# Execute run
statements.append(Forward(run("run") + Word(alphas, alphanums + "_")("run_name"))("exec_run"))
grammar_help.append(
    help_dict_create(
        name="run",
        category="Runs",
        command="run <run_name>",
        description="Execute a previously recorded run. This will execute every command and continue regardless of any failures.",
        note=NOTE_RUNS,
    )
)

# List runs
statements.append(Forward(lister + runs("runs"))("list_runs"))
grammar_help.append(
    help_dict_create(
        name="list runs",
        category="Runs",
        command="list runs",
        description="List all runs saved in the current workspace.",
        note=NOTE_RUNS,
    )
)

# Display run
statements.append(Forward(d_isplay + run("run") + Word(alphas, alphanums + "_")("run_name"))("display_run"))
grammar_help.append(
    help_dict_create(
        name="display run",
        category="Runs",
        command="display run <run_name>",
        description="Display the commands stored in a certain run.",
        note=NOTE_RUNS,
    )
)

# endregion

##########################################################################
# region - Utility
##########################################################################

# Display data
# MAJOR-RELEASE-TODO:
# In a Notebook, `x = %openad display data 'file.csv'` returns data, which
# is inconsistent with other display commands like `display mol dopamine`
statements.append(Forward(d_isplay + data("data") + desc("file_path"))("display_data"))
grammar_help.append(
    help_dict_create(
        name="display data",
        category="Utility",
        command="display data '<filename.csv>'",
        description="Display data from a csv file.",
    )
)

# --> result save --> Save data as csv
statements.append(Forward(result + save + Optional(a_s + desc("file_path")))("display_data__save"))
grammar_help.append(
    help_dict_create(
        name="save",
        category="Utility",
        command="result save [as '<filename.csv>']",
        description="Save table data to csv file.",
        parent="display data",
    )
)

# --> result open --> Explore data in browser
statements.append(
    Forward(result + CaselessKeyword("open") + Optional(CaselessKeyword("-d")("as_data")))("display_data__open")
)
grammar_help.append(
    help_dict_create(
        name="open",
        category="Utility",
        command="result open",
        description="""Explore table data in the browser.
        if you append <cmd>-d</cmd> to the end of the command <cmd>result open -d</cmd> display will result to data viewer.
        """,
        parent="display data",
    )
)

# --> result edit --> Edit data in browser
statements.append(
    Forward(result + CaselessKeyword("edit") + Optional(CaselessKeyword("-d")("as_data")))("display_data__edit")
)
grammar_help.append(
    help_dict_create(
        name="edit",
        category="Utility",
        command="result edit",
        description="""Edit table data in the browser.
        if you append <cmd>-d</cmd> to the end of the command <cmd>result open -d</cmd> display will result to data viewer.
        """,
        parent="display data",
    )
)

# --> result copy --> Copy data to clipboard, formatted for spreadheet
statements.append(Forward(result + CaselessKeyword("copy"))("display_data__copy"))
grammar_help.append(
    help_dict_create(
        name="copy",
        category="Utility",
        command="result copy",
        description="Copy table data to clipboard, formatted for spreadheet.",
        parent="display data",
    )
)

# --> result display --> Display the result in the CLI/Notebook
statements.append(
    Forward(result + CaselessKeyword("display") + Optional(CaselessKeyword("-d")("as_data")))("display_data__display")
)
grammar_help.append(
    help_dict_create(
        name="display",
        category="Utility",
        command="result display",
        description=f"""Display the result in {"Jupyter Notebook" if is_notebook_mode() else "the CLI"}.
      
        if you append <cmd>-d</cmd> to the end of the command <cmd>result open -d</cmd> display will result to data viewer.
        """,
        parent="display data",
    )
)

# --> result as dataframe --> Return the result as a dataframe
statements.append(Forward(result + a_s + CaselessKeyword("dataframe"))("display_data__as_dataframe"))
grammar_help.append(
    help_dict_create(
        name="as dataframe",
        category="Utility",
        command="result as dataframe",
        description="Return the result as dataframe (only for Jupyter Notebook)",
        parent="display data",
    )
)

# Show data
# Note: hidden from help commands because the dataviewer is not yet implemented into the GUI.
# Until then you need to make a detour via `display data` and then `result open`.
statements.append(Forward(show + data("data") + desc("file_path"))("show_data"))
# grammar_help.append(
#     help_dict_create(
#         name="show data",
#         category="Utility",
#         command="show data '<filename.csv>'",
#         description="Explore CSV data in the browser.",
#     )
# )

# Edit config file (CLI-only)
if not is_notebook_mode():
    statements.append(
        Forward(
            CaselessKeyword("edit")
            + CaselessKeyword("config")
            + desc("json_file")
            + Optional(CaselessKeyword("schema") + desc("schema"))
        )("edit_config")
    )
    grammar_help.append(
        help_dict_create(
            name="edit config",
            category="Utility",
            command="edit config '<json_config_file>' [ schema '<schema_file>']",
            description="Edit any JSON file in your workspace directly from the CLI. If a schema is specified, it will be used for validation and documentation.",
        )
    )

# endregion

##########################################################################
# region - GUI
##########################################################################

# Install gui
statements.append(Forward(install + gui)("install_gui"))
grammar_help.append(
    help_dict_create(
        name="install gui",
        category="GUI",
        command="install gui",
        description="Install the OpenAD GUI (graphical user interface).\n\nThe graphical user interface allows you to browse your workspace and visualize your datasets and molecules.",  # Partly repeated. Move to msgs()",
    )
)

# Launch gui
statements.append(Forward(launch + gui)("launch_gui"))
grammar_help.append(
    help_dict_create(
        name="launch gui",
        category="GUI",
        command="launch gui",
        description="Launch the OpenAD GUI (graphical user interface).",
    )
)

# Launch individual modules
statements.append(Forward(launch + Word(printables)("path"))("launch_gui"))
statements.append(Forward(launch + filebrowser("path"))("launch_gui"))
# grammar_help.append(
#     help_dict_create(
#         name="launch filebrowser",
#         category="GUI",
#         command="launch filebrowser",
#         description="Launch the file browser GUI module.",
#     )
# )
# grammar_help.append(
#     help_dict_create(
#         name="launch molviewer",
#         category="GUI",
#         command="launch molviewer",
#         description="Launch the molecule viewer GUI module.",
#     )
# )

# Restart gui (mostly useful for development)
statements.append(Forward(restart + gui)("restart_gui"))
grammar_help.append(
    help_dict_create(
        name="restart gui",
        category="GUI",
        command="restart gui",
        description="Terminate and then restart the GUI server.",
    )
)

# Exit gui
statements.append(Forward(q_uit + gui)("quit_gui"))
grammar_help.append(
    help_dict_create(
        name="quit gui",
        category="GUI",
        command="quit gui",
        description="Terminate the GUI server.",
    )
)


# endregion

##########################################################################
# region - LLMs
##########################################################################

# Tell me chatbot
statements.append(
    Forward(
        CaselessKeyword("tell")
        + CaselessKeyword("me")
        + ZeroOrMore(Word(alphas, alphanums + "_" + "?" + "." + " " + "," + "'" + "-" + "*" + "@" + ">"))("Chat_String")
    )("how_do_i")
)
grammar_help.append(
    help_dict_create(
        name="tell me",
        category="LLM",
        command="tell me <how to do xyz>",
        description="Ask your AI assistant how to do anything in OpenAD.",
    )
)

statements.append(
    Forward(
        CaselessKeyword("set")
        + CaselessKeyword("llm")
        + ZeroOrMore(Word(alphas, alphanums + "_" + "?" + "." + " " + "," + "'"))("llm_name")
    )("set_llm")
)
grammar_help.append(
    help_dict_create(
        name="set llm",
        category="LLM",
        command="set llm  <language_model_name>",
        description="Set the target language model name for the <cmd>tell me</cmd> command.",
    )
)

statements.append(
    Forward(CaselessKeyword("clear") + CaselessKeyword("llm") + CaselessKeyword("auth"))("clear_llm_auth")
)
grammar_help.append(
    help_dict_create(
        name="clear llm auth",
        category="LLM",
        command="clear llm auth",
        description="Clear the language model's authentication file.",
    )
)


# endregion

##########################################################################
# region - File System
##########################################################################

# List files
statements.append(
    Forward(lister + CaselessKeyword("files") + Optional(Word(alphanums + "_", alphanums + "_" + "/")("path")))(
        "list_files"
    )
)
grammar_help.append(
    help_dict_create(
        name="list files",
        category="File System",
        command="list files [ path ]",
        description="List al directories and files in your current workspace.",
    )
)

# Import file
statements.append(
    Forward(
        CaselessKeyword("import")
        + CaselessKeyword("from")
        + desc("source")
        + CaselessKeyword("to")
        + desc("destination")
    )("import_file")
)
grammar_help.append(
    help_dict_create(
        name="import",
        category="File System",
        command="import from '<external_source_file>' to '<workspace_file>'",
        description="Import a file from outside OpenAD into your current workspace.",
    )
)

# Export file
statements.append(
    Forward(
        CaselessKeyword("export")
        + CaselessKeyword("from")
        + desc("source")
        + CaselessKeyword("to")
        + desc("destination")
    )("export_file")
)
grammar_help.append(
    help_dict_create(
        name="export",
        category="File System",
        command="export from '<workspace_file>' to '<external_file>'",
        description="Export a file from your current workspace to anywhere on your hard drive.",
    )
)

# Copy file
statements.append(
    Forward(
        CaselessKeyword("copy") + CaselessKeyword("file") + desc("source") + CaselessKeyword("to") + desc("destination")
    )("copy_file")
)
grammar_help.append(
    help_dict_create(
        name="copy",
        category="File System",
        command="copy file '<workspace_file>' to '<other_workspace_name>'",
        description="Export a file from your current workspace to another workspace.",
    )
)

# Remove file
statements.append(Forward(CaselessKeyword("remove") + desc("file"))("remove_file"))
grammar_help.append(
    help_dict_create(
        name="remove",
        category="File System",
        command="remove '<filename>'",
        description="Remove a file from your current workspace.",
    )
)

# Open file
statements.append(Forward(o_pen("open") + desc("file"))("open_file"))  # From molset file
grammar_help.append(
    help_dict_create(
        name="open",
        category="File System",
        command="open '<filename>'",
        description=f"""Open a file or dataframe { 'in your browser' if is_notebook_mode() else 'in an iframe' } 

Examples:
- <cmd>open 'base_molecules.sdf'</cmd>
- <cmd>open my_dataframe</cmd>
""",
    )
)

# endregion

##########################################################################
# region - Help
##########################################################################

# Display intro.
statements.append(Forward(CaselessKeyword("intro"))("intro"))
grammar_help.append(
    help_dict_create(
        name="intro", category="Help", command="intro", description="Display an introduction to the OpenAD CLI."
    )
)
# Open documentation webpage.
statements.append(Forward(CaselessKeyword("docs"))("docs"))
grammar_help.append(
    help_dict_create(name="docs", category="Help", command="docs", description="Open the documentation webpage.")
)
# List available commands
# Note - this is controlled directly from do_help.
grammar_help.append(
    help_dict_create(name="help", category="Help", command="?", description="List all available commands.")
)

# # TO BE IMPLEMENTED LATER
# # Open advanced help.
# # Note - this is controlled directly from do_help.
# grammar_help.append(help_dict_create(
#     name="advanced help",
#     category='Help',
#     command="??",
#     description='Displays interactive help interface in the CLI.'
# ))

# Display command help
# Note - this is controlled directly from do_help.
grammar_help.append(
    help_dict_create(
        name="command help 1",
        category="Help",
        command='? ...<soft>   --> List all commands containing "..."</soft>',
        description="",
    )
)
grammar_help.append(
    help_dict_create(
        name="command help 2",
        category="Help",
        command='... ?<soft>   --> List all commands starting with "..."</soft>',
        description="",
    )
)
service_catalog_grammar(statements=statements, help=grammar_help)
"""try:
    service_catalog = get_cataloged_service_defs()
    service_grammar_add(statements=statements, help=grammar_help, service_catalog=service_catalog)
except Exception as e:
    print(e)
    pass
"""

# endregion

##########################################################################
# region - Development
# The commands in this section are not intended for general use,
# and are not documented by the help system.
##########################################################################

# Launches the demo flask app.
statements.append(Forward(CaselessKeyword("flask") + CaselessKeyword("example"))("flask_example"))

# Expose the cmd pointer
statements.append(Forward(CaselessKeyword("cmd_pointer"))("cmd_pointer"))

# endregion

# Define The Concepts of Jobs
# statements.append( Forward(lister +jobs('job') )('list_jobs'))

# This is a space that will require further thought whether to run jobs
# asynchronously from the client and how to track and save their results

# Packaging Section
orig_statements = statements.copy()
statements_def = Forward()

ii = 0
for i in statements:
    ii += 1
    statements_def |= i
statements_zom = ZeroOrMore(statements_def)


# Used to package up statements
def create_statements(cmd_pointer):
    """Create staments from Toolkits"""
    # if cmd_pointer.toolkit_current is None: # move down further to adress unset issue
    #    return
    # global statements_zom

    cmd_pointer.current_statements = orig_statements.copy()
    cmd_pointer.current_statement_defs = Forward()
    service_statements = []
    try:
        service_catalog = get_cataloged_service_defs()
        temp_help = []

        # Add model services grammar
        service_grammar_add(statements=cmd_pointer.current_statements, help=temp_help, service_catalog=service_catalog)

        # Add small molecule grammar
        smol_grammar_add(statements=cmd_pointer.current_statements, grammar_help=temp_help)

        # cmd_pointer.current_statements.extend(service_statements)

        cmd_pointer.current_help.help_model_services.clear()
        cmd_pointer.current_help.help_model_services.extend(temp_help)
        cmd_pointer.current_help.reset_help()
    except Exception as e:
        print(e)
        pass

    # Add plugin commands
    for stmt in cmd_pointer.plugins_statements:
        cmd_pointer.current_statements.append(stmt)
    cmd_pointer.current_help.help_plugins = []  # Prevent duplication
    cmd_pointer.current_help.help_plugins.extend(cmd_pointer.plugins_help)
    cmd_pointer.current_help.reset_help()

    # Add toolkit commands
    if cmd_pointer.toolkit_current is not None:
        for stmt in cmd_pointer.toolkit_current.methods_grammar:
            cmd_pointer.current_statements.append(stmt)

    # Add plugin overview commands
    for plugin_instance in cmd_pointer.plugin_instances:
        plugin_metadata = plugin_instance.metadata
        if plugin_metadata:
            plugin_namespace = plugin_metadata.get("namespace", "")
            plugin_name = plugin_metadata.get("name", "").lower()
            cmd_pointer.current_statements.append(Forward(CaselessKeyword(plugin_namespace))(plugin_namespace))
            cmd_pointer.current_statements.append(Forward(CaselessKeyword(plugin_name))(plugin_namespace))

    # Add toolkit overview commands
    for tk in _all_toolkits:
        cmd_pointer.current_statements.append(Forward(CaselessKeyword(tk))(tk))

    # Rebuild parser with updated statements
    for stmt in cmd_pointer.current_statements:
        cmd_pointer.current_statement_defs |= stmt

    # statements_zom = ZeroOrMore(statements_def)


def or_builder(options: list) -> str:
    """build or component of statement"""
    expression = "("
    the_or = ""
    for i in options:
        expression = expression + the_or + i
        the_or = "|"

    return expression + ")"


def from_builder(options: list) -> str:
    """build from clause component of statements"""
    clause = "CaselessKeyword('from')+"
    option_exp = []
    if not is_notebook_mode() and "dataframe" in options:
        options.remove("dataframe")

    for i in options:
        if i == "file":
            option_exp.append('(Suppress(CaselessKeyword("file"))+desc("from_file"))')
        elif i == "dataframe":
            option_exp.append(
                '(Suppress(CaselessKeyword("dataframe"))+' + 'Word(alphas, alphanums + "_")("from_dataframe"))'
            )
        elif i == "list":
            option_exp.append(
                '(Suppress(CaselessKeyword("list"))+ Group(Suppress("[")+delimitedList(desc)("from_list")+Suppress("]")))'
            )
    if len(option_exp) > 0:
        clause = clause + or_builder(option_exp)

    # elif len(option_exp) ==1:
    #    clause=clause+option_exp[0]
    else:
        raise ValueError("invalid 'From' Clause Structure")

    return clause


def statement_builder(toolkit_pointer, inp_statement):
    """builds statements from toolkit function defintions"""
    #####################################################################
    #
    # This section deals with Method call like statements, they will always be prefixed with 'exec'.
    ########################################################################################################
    if inp_statement["exec_type"] == "method":
        expression = "Suppress(e_xec)+" + 'CaselessKeyword ("' + inp_statement["command"] + '")'
        expression = (
            expression
            + '+ Suppress("(") +'
            + optional_parameter_list(inp_statement, "fixed_parameters")
            + ' +Suppress(")"))'
        )

    elif inp_statement["exec_type"] == "standard_statement":
        key_leading_words = inp_statement["command"].split()
        expression = ""
        plus = ""

        for i in key_leading_words:
            expression = expression + plus + 'CaselessKeyword ("' + i + '")'
            plus = "+"

        if "SINGLE_PARM" in inp_statement:
            if len(inp_statement["SINGLE_PARM"]) > 0:
                expression = expression + "+" + actual_parameter_list(inp_statement, "SINGLE_PARM")

        if "from" in inp_statement:
            if len(inp_statement["from"]) > 0:
                expression = expression + "+" + (from_builder(inp_statement["from"])) + "('from_source')"
        if "USING" in inp_statement:
            if inp_statement["USING"] is not None and len(inp_statement["USING"]) > 0:
                expression = (
                    expression
                    + '+ Optional( (CaselessKeyword ("USING")+ Suppress("(") +'
                    + optional_parameter_list(inp_statement, "USING")
                    + '+Suppress(")") )("USING"))'
                )
        if "RETURN_AS_DATA" in inp_statement:
            expression = (
                expression
                + '+ Optional(Suppress(CaselessKeyword ("return")) + Suppress(CaselessKeyword ("as"))+Suppress(CaselessKeyword ("data")))("return_as_data") '
            )

        if "SAVE_AS" in inp_statement:
            expression = (
                expression
                + '+ Optional(Suppress(CaselessKeyword ("save")) + Suppress(CaselessKeyword ("as"))+desc("results_file") )("save_as") '
            )
        if "use_saved" in inp_statement and len(inp_statement["use_saved"]) > 0:
            if inp_statement["use_saved"] == "True":
                expression = expression + "+" + "Optional(CaselessKeyword('use_saved'))('use_saved')"
        if "async" in inp_statement and len(inp_statement["async"]) > 0:
            if inp_statement["async"] == "both":
                expression = expression + "+" + "Optional(CaselessKeyword('async'))('do_async')"
            if inp_statement["async"] == "only":
                expression = expression + "+" + "CaselessKeyword('async')"

        expression = expression + ")"

    #####################################################################
    #
    # This section deals with search type statements.
    ########################################################################################################
    elif inp_statement["exec_type"] == "search_statement":
        expression = ""
        a_char = " "
        for i in inp_statement["command"].split():
            if i == "<in_clause>" and "in_clause" in inp_statement:
                for in_clause in inp_statement["in_clause"]:
                    expression = expression + a_char + ' CaselessKeyword ("' + in_clause + '")'
                    if inp_statement["in_clause"][in_clause] == "desc":
                        expression = expression + a_char + ' desc("' + in_clause + '")'
                    if inp_statement["in_clause"][in_clause] == "str":
                        expression = expression + a_char + ' Word(alphas, alphanums + "_")(' + in_clause + ")"
            else:
                expression = expression + a_char + ' CaselessKeyword ("' + i + '")'

            if a_char == " ":
                a_char = " + "

        expression = expression + " +desc('val')('" + inp_statement["subject"] + "') "

        ####################################################################################################
        #
        # "USING" clase specifies any indexes, attributes or other services that optionally can be used
        # currently if a field is mandatory this will need to be picked up in the toolkit program
        # e.g. Using index_type="arxiv" contstraint=chemicals maxhops=1
        ####################################################################################################
        if "USING" in inp_statement:
            if inp_statement["USING"] is not None and len(inp_statement["USING"]) > 0:
                expression = (
                    expression
                    + '+  Optional(CaselessKeyword ("USING")+ Suppress("(") +'
                    + optional_parameter_list(inp_statement, "USING")
                    + '+Suppress(")") )("USING")'
                )
        ######################################################################################################
        #
        # Show (<data_types>) is to specify what data to retrieve.
        # This is optional clause, any mandatory values need to be managed by the toolkit.
        # ####################################################################################################

        if "SHOW" in inp_statement:
            if inp_statement["SHOW"] is not None:
                options = []
                for i in inp_statement["SHOW"]:
                    options.append(i)

                expression = (
                    expression
                    + '+ Optional(CaselessKeyword ("SHOW") +Suppress("(")+'
                    + oneormore_str(options)
                    + ' +Suppress(")"))("show_data") '
                )
        ######################################################################################################
        #
        # Estimate Only is for statements that may take a long time and there is an estimate capability.
        # ####################################################################################################
        if "ESTIMATE_ONLY" in inp_statement:
            expression = (
                expression
                + '+ Optional(Suppress(CaselessKeyword ("ESTIMATE")) +Suppress(CaselessKeyword ("ONLY")) )("estimate_only") '
            )
        if "RETURN_AS_DATA" in inp_statement:
            expression = (
                expression
                + '+ Optional(Suppress(CaselessKeyword ("return")) + Suppress(CaselessKeyword ("as"))+Suppress(CaselessKeyword ("data")))("return_as_data") '
            )

        if "SAVE_AS" in inp_statement:
            expression = (
                expression
                + '+ Optional(Suppress(CaselessKeyword ("save")) + Suppress(CaselessKeyword ("as"))+desc("results_file") )("save_as") '
            )

        expression = expression + ")"

    try:
        # There is no alternative for eval for this logic
        toolkit_pointer.methods_grammar.append(
            eval(" Forward( " + expression + ' ("toolkit_exec_' + inp_statement["command"] + '")')
        )

        toolkit_pointer.methods.append(inp_statement["command"])

        toolkit_pointer.methods_execute.append(inp_statement["method"])

        toolkit_pointer.methods_library.append(inp_statement["library"])

        toolkit_pointer.methods_dict.append(inp_statement)

        # TEMPORARY toolkit support
        # We switched help_dict_create to store command string and aliases in "commands" field instead of "command".
        # We need to translate the toolkit help to the new format.
        #
        # Before: toolkit_pointer.methods_help.append(inp_statement["help"])
        help_statement = inp_statement["help"]
        cmd = help_statement.pop("command")
        help_statement["commands"] = [cmd]
        #
        toolkit_pointer.methods_help.append(help_statement)

    except Exception as err:
        fwd_expr = "Forward( " + expression + ' ("toolkit_exec_' + inp_statement["command"] + '")'
        output_error(msg("err_add_command", inp_statement["command"], fwd_expr, err))

    return True


def oneormore_str(options: list):
    """create one or more parsing string"""
    expression = " OneOrMore("
    divider = " "

    for i in options:
        expression = expression + divider + ' CaselessKeyword ("' + i.strip() + '") '

        divider = "|"

    expression = expression + ")"
    return expression


def optional_parameter_list(inp_statement: dict, clause: str):
    """Create an optional parameter list for a clause"""
    ii = 0
    expression = " "
    for i in inp_statement[clause]:
        if ii == 0:
            expression = expression + " "
            if inp_statement[clause][i] == "str":
                expression = (
                    expression
                    + " ZeroOrMore(Group( CaselessKeyword ('"
                    + i
                    + "') +Suppress('=')+key_val_expr('val'))('"
                    + i
                    + "'))"
                    + " "
                )
            elif inp_statement[clause][i] == "desc":
                expression = (
                    expression
                    + " ZeroOrMore(Group( CaselessKeyword ('"
                    + i
                    + "') +Suppress('=')+desc('val'))('"
                    + i
                    + "'))"
                    + " "
                )
            else:
                expression = (
                    expression
                    + " ZeroOrMore(Group( CaselessKeyword ('"
                    + i
                    + "') +Suppress('=')+key_val_expr_num('val'))('"
                    + i
                    + "'))"
                    + " "
                )
        else:
            if inp_statement[clause][i] == "str":
                expression = (
                    expression
                    + " & ZeroOrMore(Group( CaselessKeyword ('"
                    + i
                    + "') +Suppress('=')+key_val_expr('val'))('"
                    + i
                    + "'))"
                    + " "
                )
            elif inp_statement[clause][i] == "desc":
                expression = (
                    expression
                    + " & ZeroOrMore(Group( CaselessKeyword ('"
                    + i
                    + "') +Suppress('=')+desc('val'))('"
                    + i
                    + "'))"
                    + " "
                )
            else:
                expression = (
                    expression
                    + " & ZeroOrMore(Group( CaselessKeyword ('"
                    + i
                    + "') +Suppress('=')+key_val_expr_num('val'))('"
                    + i
                    + "'))"
                    + " "
                )
        ii = 1

    return expression


def actual_parameter_list(inp_statement: dict, clause: str):
    """create parameter list for clause"""
    ii = 0
    expression = " "
    for i in inp_statement[clause]:
        if ii == 0:
            expression = expression + " "
            if inp_statement[clause][i] == "str":
                expression = expression + " key_val_expr('" + i + "') "
            elif inp_statement[clause][i] == "desc":
                expression = expression + "desc('" + i + "') "
            else:
                expression = expression + " key_val_expr_num('" + i + "') "
        else:
            if inp_statement[clause][i] == "str":
                expression = expression + ", key_val_expr('" + i + "') "
            elif inp_statement[clause][i] == "desc":
                expression = expression + ",desc('" + i + "') "
            else:
                expression = expression + ", key_val_expr_num('" + i + "') "
        ii = 1
    return expression


def output_train_statements(cmd_pointer):
    """Create training statements for the Tell Me Prompt"""
    training_statements = []
    i = 0
    try:
        if not os.path.exists(os.path.expanduser(os.path.expanduser(cmd_pointer.home_dir + "/prompt_train/"))):
            os.mkdir(os.path.expanduser(os.path.expanduser(cmd_pointer.home_dir + "/prompt_train/")))
    except Exception as err:
        raise FileExistsError(
            f"unable to create prompt soc directory {str(os.path.expanduser(cmd_pointer.home_dir+'/prompt_train/'))}::: {err}"
        ) from err

    for training_file in glob.glob(
        os.path.expanduser(str(os.path.expanduser(cmd_pointer.home_dir + "/prompt_train/")) + "/*")
    ):
        os.remove(training_file)

    while i < len(grammar_help):
        training_statements.append(
            {
                "command_group": "base",
                "command_name": grammar_help[i]["name"].replace("_", " "),
                "command_syntax": tags_to_markdown("\n" + "\n".join(grammar_help[i]["commands"])),
                "command_help": _parse_description(grammar_help[i]["description"]),
            }
        )
        i += 1
    i = 0
    while i < len(cmd_pointer.current_help.help_model_services):
        training_statements.append(
            {
                "command_group": "base",
                "command_name": cmd_pointer.current_help.help_model_services[i]["name"].replace("_", " "),
                "command_syntax": tags_to_markdown(
                    "\n" + "\n".join(cmd_pointer.current_help.help_model_services[i]["commands"]),
                ),
                "command_help": _parse_description(
                    cmd_pointer.current_help.help_model_services[i]["description"],
                ),
            }
        )
        i += 1

    training_file = open(
        os.path.expanduser(cmd_pointer.home_dir + "/prompt_train/base_commands.cdoc"),
        "w",
        newline="\n",
        encoding="utf-8",
    )
    training_file.write(
        """openad client information

        For information in providing answers to how to or Help questions from users :

        The below describes openad clients domain specific language (DSL) for managing science activities using the DSL

        Vocabulary:
            DSL: Domain Specific Language or DSL  that is implemented for the openad client and all commands are formatted in
            Set: set the current status of a system variable e.g. Workspace or Context(current toolkit )
            Login: login into a Plugins target system if requires manual login
            Logout: logout of a system if required
            List: list a set of objects
            Get: get the details of an object
            Create: create an object (e.g. workspace or run)
            Search: a repository or object
            Exec/Execute: execute a function
            Display: display molecule, file or result set
            Show:  Show a data set using a utility that enables you to manipulate or  diagrammatically view it.
            Backup: backup a plugin or workspace
            Add: add a function or plugin
            Remove: delete an object
            Save: Save a run or file of some kind
            Load: load a file from project directory to Target system
            pyparsing_statement: a statement defined using pyparsing for the domain specific language
            help_text: description of the Domain Specific language statement defined in a pyparsing_statement
            toolkit: these are contextual plugins that are available one at a time for providing specific functionality to the user. Valid toolkits are DS4SD (deep Search),  RXN (retro synthesis), ST4SD(simulation toolkit)
                   The Deep Search toolkit and RX toolkits have separate help outlining specific commands avilable to the user
            History: History of DSL commands for a given Workspace
            run: list of sequential commands saved by the user')
            molecule working set: is a set of molecules in memory that can added to using the 'add molecule' command  and also loaded from a molecule-set and maipulated by commands suchs as 'display molecule', 'add Molecule','create molecule', 'remove molecule' 'merge mol-set'
            molecules in the molecule working set contain the following data
                - identifiers shuch as names, synonyms, inchi , inchikey ,canonical smiles, isomeric smiles and the CID or compound ID sourced from pubchem as an identfier
                - properties : data points such as atom_stereo_count, bond_stereo_count, canonical_smiles, charge, cid, complexity, conformer_count_3d, conformer_id_3d, conformer_model_rmsd_3d, conformer_rmsd_3d, coordinate_type, covalent_unit_count, defined_atom_stereo_count, defined_bond_stereo_count, effective_rotor_count_3d, exact_mass, feature_acceptor_count_3d, feature_anion_count_3d, feature_cation_count_3d, feature_count_3d, feature_donor_count_3d, feature_hydrophobe_count_3d, feature_ring_count_3d, h_bond_acceptor_count, h_bond_donor_count, heavy_atom_count, inchi, inchikey, isomeric_smiles, isotope_atom_count, iupac_name, mmff94_energy_3d, mmff94_partial_charges_3d, molecular_formula, molecular_weight, monoisotopic_mass, multipoles_3d, multipoles_3d, pharmacophore_features_3d, pharmacophore_features_3d, rotatable_bond_count, tpsa, undefined_atom_stereo_count, undefined_bond_stereo_count, volume_3d, x_steric_quadrupole_3d, xlogp, y_steric_quadrupole_3d, z_steric_quadrupole_3d
                - synonyms : names that it is known by both commerically and scientifically
                - analysis : records of selected result sets from supporting analysis functions such as predict retrosynthesis
            molecules when displayed will show these data points.
            molecule-set: a list of molecule dictionaries that has been stored on disk under a molecule-set name and can be loaded into the molecule working set in a users sessions
            The short form of 'molecule-set' is 'molset'
            The short form of 'molecule' is 'mol' 
            
            The Model Service is a capability to register and launch model services for property prediction and data set generation and allows you to launch ones you catalog yourself or remotely catalog already running services.

           \@
            If a user asks for parameters or options this refers to the parameters that can be given to a function. Make sure all parameters are provided to the user
            
            The Following commands are used to work with the molecule working set of molecules:
                - add molecule <name> | <smiles> | <inchi> | <inchikey> | <cid>   [as '<name>' ] [ basic ] [ force ]
                - display molecule <name> | <smiles> | <inchi> | <inchikey> | <cid>
                - rename molecule <molecule_identifer_string> as <molecule_name>
                - export molecule <name> | <smiles> | <inchi> | <inchikey> | <cid> [ as file ]
                - remove molecule <name> | <smiles> | <inchi> | <inchikey> | <cid>
                - list molecules
                - clear molecules
                - load molecules using file '<csv_or_sdf_filename>' [ merge with pubchem ]
                - load molecules using dataframe <dataframe>  [merge with pubchem]
                - export molecules [ as <csv_filename> ]
                - show molecule <name> | <smiles> | <inchi> | <inchikey> | <cid>
                - show molecules using ( file '<mols_file>' | dataframe <dataframe> ) [ save as '<sdf_or_csv_file>' | as molsobject ]
                - save molecule-set as <molecule_set_name>
                - load molecule-set|molset <molset_name>
                - merge molecule-set|molset <molset_name> [merge only] [append only]
                - list molecule-sets
                - enrich molecules with analysis
                - @(<name> | <smiles> | <inchi> | <inchikey> | <cid>)>><molecule_property_name>
            
            A parameter designated by <smiles> refers to a parameter with a caluse such as C(C(C1C(=C(C(=O)O1)O)O)O)O 

            ?: will display help and if positioned prior to a command will display help options for that command
        \@
             The Model Service is a capability to register and launch model services for property generation and data set generation and allows you to launch ones you catalog yourself or remotely catalog already running services.
    
            ***Molecules Properties***

        To generate a molecules Property we use the 'GET MOLECULE PROPERTY" command 
        
            ***GET MOLECULE PROPERTY COMMAND***
            `GET MOLECULE PROPERTY` is a standard command that is used for different model commands. when a user catalogs a model that generates properties 1 or commands are created starting with the namespace or prefix for the model e.g. prop and the syntax for the different available commands e.g. for a 
            
            Question: What is the command Syntax for generating a molecules property ?
            Answer: {
            The format 'GET MOLECULE PROPERTY" command is as follows: Command Syntax: <cmd><model_prefix> GET MOLECULE PROPERTY @mols | <molecule property identifer> | [ <list of molecule_property identifers> ] FOR  <smiles_string> | [<list of smiles strings >] USING(<option>=<value>) (merge with molecules|mols) </cmd>
           
            - `<model_prefix>` : Nomimal user defined name provided to a model service when it is cataloged
            - `<molecule_property_identifer>` : valid property idenfitier for the given model service that is to be generated in the request like 'esol','xlogp' AA1R
            - `[<list of molecule property identifer>]` : valid list of property idenfitiers in square brackets and comma separated e.g. [ qes, esol ] for the given model service that is to be generated in the request
            - `@mols`  specifies the list in the current molecules working set.
            - `<smiles string>` : a string denoting a valid canonical smiles like 'C(C(C1C(=C(C(=O)O1)O)O)O)O' that is either free text or in single quotes.
            - `[ <list of smiles strings> ]` :  a list of strings denoting a valid canonical smiles that is either free text or in single quotes and separated by commas enclided in square brackets
            The clause `merge with mols` will merge the resulting molecule properties with the memory molecule working set.
            The Using clause is incased in normal brackets () and composed of a space delimited string of valid <options>  and their values  e.g. <option>=<value>

            Lists of smiles or properties are contained in square brackets, smiles strings can be in single quoted strings or without quotes, hoever if special characters are in the strings it is best to place them in quotes.
            The Using clause is optional and properties are grouped into commands bsed on common properties.
            
            With a model_prefix called `prop` the following molecule statements are examples:
            
            generate the esol property FOR a list of molecules defined by smiles string
            - <cmd>prop GET MOLECULE PROPERTY esol FOR ['C(C(C1C(=C(C(=O)O1)O)O)O)O','[H-]']</cmd>
            - <cmd>prop GET MOLECULE PROPERTY qed FOR ['C(C(C1C(=C(C(=O)O1)O)O)O)O','[H-]']</cmd>
             - <cmd>prop GET MOLECULE PROPERTY lipinski FOR ['C(C(C1C(=C(C(=O)O1)O)O)O)O','[H-]']</cmd>

            generate a list of properties FOR a single smiles string
            - <cmd>prop GET MOLECULE PROPERTY [qed,esol] FOR 'C(C(C1C(=C(C(=O)O1)O)O)O)O'</cmd>

            Generate a list of properties FOR a list of smiles strings
            - <cmd>prop GET MOLECULE PROPERTY [qed,esol] FOR [ C(C(C1C(=C(C(=O)O1)O)O)O)O ,[H-] ]</cmd>

            Generate a list of properties FOR the smiles strings in the molecule working set in memory
            - <cmd>prop GET MOLECULE PROPERTY [qed,esol] FOR @mols </cmd>

            Generate a list of properties FOR the smiles strings in the molecule working set in memory and merge it back into the molecule working set
            - <cmd>prop GET MOLECULE PROPERTY [qed,esol] FOR @mols merge with mols</cmd>

            Generate a single property for a single smiles string
            - <cmd>prop GET MOLECULE PROPERTY esol FOR C(C(C1C(=C(C(=O)O1)O)O)O)O</cmd>
            - <cmd>prop GET MOLECULE PROPERTY qed FOR [H-] </cmd>
            - <cmd>prop GET MOLECULE PROPERTY xlogp FOR [H-] </cmd>
            - <cmd>prop GET MOLECULE PROPERTY lipinski FOR [H-] </cmd>

            Generate a single property for a single smiles string and merge it back into the molecule working set
            - <cmd>prop GET MOLECULE PROPERTY esol FOR C(C(C1C(=C(C(=O)O1)O)O)O)O merge with mols</cmd> 
            - <cmd>prop GET MOLECULE PROPERTY qed FOR [H-] merge with mols</cmd> 
            - <cmd>prop GET MOLECULE PROPERTY xlogp FOR [H-] merge with mols</cmd> 
            - <cmd>prop GET MOLECULE PROPERTY lipinski FOR [H-] merge with mols</cmd> 

            Generate A SINGLE molecules property FOR a list of smiles strings using an additional parameter
            - <cmd>prop GET MOLECULE PROPERTY activity_against_target FOR ['C(C(C1C(=C(C(=O)O1)O)O)O)O','[H-]'] using(target=drd2)</cmd>

            Generate A SINGLE molecules property FOR a list of smiles strings using an additional parameter and merge it back into the molecule working set
            - <cmd>prop GET MOLECULE PROPERTY activity_against_target FOR ['C(C(C1C(=C(C(=O)O1)O)O)O)O','[H-]'] using(target=drd2) merge with mols</cmd>
            }
            \@
            ***GET PROTEIN PROPERTY***

            when identifying a protein we can use what is known as a FASTA string. In bioinformatics and biochemistry, the FASTA format is a text-based format for representing either nucleotide sequences or amino acid (protein) sequences, in which nucleotides or amino acids are represented using single-letter codes.
            An example of a FASTA string is `MADQLTEEQIAEFKEAFSLFDKDGDGTITTKELGTVMRSLGQNPTEAELQDMINEVDADGNGTIDFPEFLTMMARKMKDTDSEEEIREAFRVFDKDGNGYISAAELRHVMTNLGEKLTDEEVDEMIREADIDGDGQVNYEEFVQMMTAK`
            - `<model_prefix>` : Nomimal user defined name provided to a model service when it is cataloged
            - `<protein_property_identifer>` : valid proitein property idenfitier such as length, boman_index, aliphaticity, hydrophobicity, aromaticity, instability
            - `[<list_of_protein_property_identifers>]` : valid list of property idenfitiers in square brackets and comma separated e.g. [boman_index, aliphaticity, hydrophobicity, aromaticity, instability]
            - `<FASTA_string>` : a string denoting a valid FASTA string  'MDITIHNPLIRRPLFSWLAPSRIFDQIFGEHLQESELLPASPSLSPFLMRPIFRMPSWLETGLSEMRLEKDKFSVNLDVKHFSPEELKVKVLGDMVEIHGKHEERQDEHGFIAREFNRKYRIPADVDPLTITSSLSLDGVLTVSAPRKQSDVPERSIPITREEKPAIAGAQRK' that is e in single quotes.
            - `[ <list of FASTA strings> ]` :  a list of strings denoting a valid FASTA representation of proteins that is in single quotes and separated by commas enclided in square brackets []
            The Using clause is incased in normal brackets () and composed of a space delimited string of valid <options>  and their values  e.g. <option>=<value>
            this is what we use for identifying proteins for property generation.
            The Command syntax for 'GET PROTEIN PROPERTY' command is as follows :<cmd><model_prefix> GET MOLECULE PROTEIN <protein property identifer> | [<list of protein property identifers>] FOR  <FASTA string> | [<list of FASTA strings> ] USING(<option>=<value>)</cmd>
            <cmd>prop get protein property [ charge_density, charge ]  for ['MKYNNRKLSFNPTTVSIAGTLLTVFFLTRLVLSFFSISLFQLVTFQGIFKPYVPDFKNTPSVEFYDLRNYQGNKDGWQQGDRILFCVPLRDASEHLPMFFNHLNTMTYPHNLIDLSFLVSDSSDNTMGVLLSNLQMAQSQQDKSKRFGNIEIYEKDFGQIIGQSFSDRHGFGAQGPRRKLMARARNWLGSVALKPYHSWVYWRDVDVETIPTTIMEDLMHHDKDVIVPNVWRPLPDWLGNIQPYDLNSWKESEGGLQLADSLDEDAVIVEGYPEYATWRPHLAYMRDPNGNPEDEMELDGIGGVSILAKAKVFRTGSHFPAFSFEKHAETEAFGRLSRRMNYNVIGLPHYVIWHIYEPSSDDLKHMAWMAEEEKRKLEEERIREFYNKIWEIGFEDVRDQWNEERDSILKNIDSTLNNKVTVDWSEEGDGSELVDSKGDFVSPNNQQQQQQQQQQQQQQQQQQQQQQLDGNPQGKPLDDNDKNKKKHPKEVPLDFDPDRN','MQYLNFPRMPNIMMFLEVAILCLWVVADASASSAKFGSTTPASAQQSDVELEPINGTLNYRLYAKKGRDDKPWFDGLDSRHIQCVRRARCYPTSNATNTCFGSKLPYELSSLDLTDFHTEKELNDKLNDYYALKHVPKCWAAIQPFLCAVFKPKCEKINGEDMVYLPSYEMCRITMEPCRILYNTTFFPKFLRCNETLFPTKCTNGARGMKFNGTGQCLSPLVPTDTSASYYPGIEGCGVRCKDPLYTDDEHRQIHKLIGWAGSICLLSNLFVVSTFFIDWKNANKYPAVIVFYINLCFLIACVGWLLQFTSGSREDIVCRKDGTLRHSEPTAGENLSCIVIFVLVYYFLTAGMVWFVFLTYAWHWRAMGHVQDRIDKKGSYFHLVAWSLPLVLTITTMAFSEVDGNSIVGICFVGYINHSMRAGLLLGPLCGVILIGGYFITRGMVMLFGLKHFANDIKSTSASNKIHLIIMRMGVCALLTLVFILVAIACHVTEFRHADEWAQSFRQFIICKISSVFEEKSSCRIENRPSVGVLQLHLLCLFSSGIVMSTWCWTPSSIETWKRYIRKKCGKEVVEEVKMPKHKVIAQTWAKRKDFEDKGRLSITLYNTHTDPVGLNFDVNDLNSSETNDISSTWAAYLPQCVKRRMALTGAATGNSSSHGPRKNSLDSEISVSVRHVSVESRRNSVDSQVSVKIAEMKTKVASRSRGKHGGSSSNRRTQRRRDYIAAATGKSSRRRESSTSVESQVIALKKTTYPNASHKVGVFAHHSSKKQHNYTSSMKRRTANAGLDPSILNEFLQKNGDFIFPFLQNQDMSSSSEEDNSRASQKIQDLNVVVKQQEISEDDHDGIKIEELPNSKQVALENFLKNIKKSNESNSNRHSRNSARSQSKKSQKRHLKNPAADLDFRKDCVKYRSNDSLSCSSEELDVALDVGSLLNSSFSGISMGKPHSRNSKTSCDVGIQANPFELVPSYGEDELQQAMRLLNAASRQRTEAANEDFGGTELQGLLGHSHRHQREPTFMSESDKLKMLLLPSK']</cmd>
 \@ \n"""
    )

    """training_file.write(
        "The following dictionaries are the Domain Specific Language (DSL) commands available in the base openad client.\n \
              The users Domain Specific Language commands will be interpreted from these command definitions which the application\
                  trainslates into the Domain specific Language (DSL) .. for reference \n command_name : is the name of the command\n command_group :\
                      is the toolkit group the command belongs to, command_syntax : is the Syntax description for the command \n command_help :\
                          is the associated help instructions and examples on how to use the command.  \@\n"
    )"""
    # for i in training_statements:
    #    training_file.write(str(i) + "\\@\n")
    training_file.close()
    cmds = []
    cmds.extend(_compile_section(organize_commands(grammar_help)))
    cmds.extend(_compile_section(organize_commands(cmd_pointer.current_help.help_model_services)))
    commands = "\n".join(cmds)
    i = 0
    new_line_replace = """
"""
    for command in cmds:
        i = i + 1
        training_file = open(
            os.path.expanduser(cmd_pointer.home_dir + f"/prompt_train/individual_command_{str(i)}.cdoc"),
            "w",
            newline=new_line_replace,
            encoding="utf-8",
        )
        command.replace("\n", new_line_replace)
        training_file.write(command)
        training_file.close()
    """import json

    training_file = open(
        os.path.expanduser(cmd_pointer.home_dir + f"/prompt_train/individual_command.cdoc"),
        "w",
        newline="\n",
        encoding="utf-8",
    )
    for i in training_statements:

        training_file.write(
            json.dumps(
                i,
                ensure_ascii=False,
                indent=4,
            )
        )"""

    for i in cmd_pointer.settings["toolkits"]:
        token, a_toolkit = toolkit_main.load_toolkit(i, for_training=True)
        training_statements = []
        training_file = open(
            os.path.expanduser(cmd_pointer.home_dir + "/prompt_train/toolkit_" + i + ".cdoc"),
            "w",
            newline="\n",
            encoding="utf-8",
        )
        training_file.write("This file represents the description of the " + a_toolkit.toolkit_name + ". \n\n")

        training_file.write(
            "this is the description of the overall capability of the " + a_toolkit.toolkit_name + ". \n \n"
        )
        training_file.write(
            "The following are the  command definitions for the "
            + a_toolkit.toolkit_name
            + " toolkit with their  command definition and description \n "
            + str(a_toolkit.toolkit_description.replace("\n", new_line_replace))
            + " \n When displaying an answer always interpret the command specified \@ \n"
        )

        x = 0
        while x < len(a_toolkit.methods):
            training_statements.append(
                {
                    "toolkit group": a_toolkit.toolkit_name,
                    "command_name": a_toolkit.methods_help[x]["name"],
                    "commands": "\n".join(a_toolkit.methods_help[x]["commands"]),
                    "command_help": a_toolkit.methods_help[x]["description"],
                }
            )
            x += 1
        for i in training_statements:
            training_file.write("\@ " + str(i))
        training_file.close()
    return


def _parse_description(description):
    # description = tags_to_markdown(description)

    # Style notes as blockquotes, and ensure they're always
    # followed by an empty line, to avoid the next lines to
    # be treated as part of the blockquote.
    description = re.sub(
        r"(\*\*Note:\*\*.+?)(\n{1,})",
        lambda match: (
            f"  > {match.group(1)}\n\n" if len(match.group(2)) == 1 else f"  > {match.group(1)}{match.group(2)}"
        ),
        description,
        flags=re.MULTILINE,
    )

    # description = description.splitlines()
    # description = "\n".join([line.strip() for line in description])
    return description.strip()


# Compile all commands of a single section.
def _compile_section(cmds_organized):
    output = []

    for category in cmds_organized:
        category_output = ""
        for cmd_str, cmd_description in cmds_organized[category]:
            category_output = (
                category_output
                + f"Command Category: {category}\n Command Syntax: `{cmd_str.strip()}`\nCommand Description: {_parse_description(cmd_description)}<br>\n \@"
            )

        output.append(category_output)
    return output


# Need to fix later
#!/usr/local/opt/python@3.9/bin/python3.9
