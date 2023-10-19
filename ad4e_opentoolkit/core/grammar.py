# Note: this file is organized in regions, which collapse in VS Code:
# - Collapse all regions:       hold cmd, then hit K followed by 1
# - Expand all regions:         hold cmd, then hit K followed by J
# - Collapse to any level:      first expand everything, then hold cmd, then hit K followed by any number


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
    # Literal,
    # replaceWith,
    # Combine,
    # pyparsing_test,
    # ParseException,
)

# Core
from ad4e_opentoolkit.core.help import help_dict_create

# Helpers
from ad4e_opentoolkit.helpers.general import is_notebook_mode
from ad4e_opentoolkit.helpers.output import output_error, msg

# Global variables
from ad4e_opentoolkit.app.global_var_lib import _all_toolkits


get, lister, description, using, create, s_et, unset, workspace, workspaces, context, jobs, e_xec, a_s, optimize, w_ith, toolkits, toolkit, GPU, experiment, add, run, save, runs, show, molecules, file, d_isplay, history, data, remove, result = \
    map(CaselessKeyword, "get list description using create set unset workspace workspaces context jobs exec as optimize with toolkits toolkit gpu experiment add run save runs show molecules file display history data remove result".split())
string_value = alphanums

##########################################################################

skobki = "(" + ZeroOrMore(CharsNotIn(")")) + ")"
field_def = OneOrMore(
    Word(alphas, alphanums + "_\"':-") | skobki)


def field_act(s, loc, tok):
    return ("".join(tok)).replace('"', '\\"')
    # return ("<" + tok[0] + "> " + " ".join(tok)).replace('"', '\\"')


field_def.setParseAction(field_act)
field_list_def = delimitedList(field_def)


def field_list_act(toks):
    return " , ".join(toks)


field_list_def.setParseAction(field_list_act)
# desc=QuotedString("'" ,escQuote='"')
desc = QuotedString("'", escQuote='\\')

# 3
# This section needs revisiting and clean up

# interpret inline 'string' as Suppress('string'),
# instead of LBRACE,EQ,RBRACE,HASH = map(Suppress, '{=}#')
ParserElement.inlineLiteralsUsing(Suppress)

# be sure to exclude special characters when using printables
# name_expr = Word(printables, excludeChars='{}')
# key_val_expr = Word(printables)
name_expr = Word(alphanums + "_" + ".")
key_val_expr = Word(alphanums + "_" + ".")
key_val_expr_num = Word(nums)
key_val_expr_alpha = Word(alphanums + "_" + ".")
# p1('name') is equivalent to p1.setResultsName('name')
# p1 | p2 is equivalent to MatchFirst(p1, p2)
# if lineEnd() matches first, there is no value.
# then use a parse action to return the string 'NONE' as value instead
# else, match a regular key_value
# also, you have to use Group because key_val_line is a repeating element
# key_val_line = Group(name_expr('key') + (lineEnd().setParseAction(lambda t 'NONE') | key_val_expr)('val'))
key_val_line = Group(name_expr('key') + Suppress('=') + key_val_expr('val'))
key_val_lines = OneOrMore(key_val_line)
complex_obj = Forward()
complex_obj << Group(name_expr('obj_name') + '{' + ZeroOrMore(key_val_lines | complex_obj) + '}')('objects')
GPU_CLAUSE = (Suppress(optimize + w_ith) + GPU)
AS_EXPERIMENT = (Suppress(a_s + experiment) + Word(alphas, alphanums + "_" + "."))
obj = Forward()
obj << Suppress("(") + ZeroOrMore(key_val_lines) + Suppress(")")


delim_value = Group(delimitedList(string_value))("list")

##########################################################################

statement = Forward()

statements = []  # Statement definitions
grammar_help = []  # Help text


##########################################################################
# region - General
##########################################################################

# General splash screen
statements.append(Forward(CaselessKeyword('openad'))('welcome'))
grammar_help.append(help_dict_create(
    name="welcome",
    category='General',
    command="openad",
    description='Display the openad splash screen.'
))

# Get status
statements.append(Forward(get + CaselessKeyword('status'))('get_status'))
grammar_help.append(help_dict_create(
    name="get status",
    category='General',
    command="get status",
    description='Display the currently selected workspace and toolkit.'
))

# Display history
statements.append(Forward(d_isplay + history('history'))('display_history'))
grammar_help.append(help_dict_create(
    name="Display History",
    category='General',
    command="display history ",
    description='Display the last 30 commands run in your current workspace.'
))

# Clear sessions
statements.append(Forward(CaselessKeyword('clear') + CaselessKeyword('sessions'))('clear_sessions'))
grammar_help.append(help_dict_create(
    name="Clear Sessions",
    category='General',
    command="clear sessions",
    description='Clear any other sessions that may be running.'
))

# endregion

##########################################################################
# region - Workspaces
##########################################################################

info_workspaces = '\n<soft>To learn more about workspaces, run <cmd>workspace ?</cmd></soft>'

# Set workspaces
statements.append(Forward(s_et + workspace('workspace') + Word(alphas, alphanums + "_")('Workspace_Name'))('set_workspace_statement'))
grammar_help.append(help_dict_create(
    name="set workspace",
    category='Workspaces',
    command="set workspace <workspace_name>",
    description=f'Change the current workspace.{info_workspaces}'
))

# Get workspace
statements.append(Forward(get + workspace('workspace') + Optional(Word(alphas, alphanums + "_")('Workspace_Name')))('get_workspace'))
grammar_help.append(help_dict_create(
    name="get workspace",
    category='Workspaces',
    command="get workspace [ <workspace_name> ]",
    description=f'Display details a workspace. When no workspace name is passed, details of your current workspace are displayed.{info_workspaces}'
))

# Create workspaces
statements.append(Forward(
    create + workspace('workspace') + Word(alphas, alphanums + "_")('Workspace_Name') +
    Optional(description('description') + Suppress("(") + desc("proj_desc") + Suppress(")"))('description_option') +
    Optional(CaselessKeyword('on') + CaselessKeyword('path') + desc("w_path"))('workspace_path')
)('create_workspace_statement'))
grammar_help.append(help_dict_create(
    name="create workspace",
    category='Workspaces',
    command="create workspace <workspace_name> [ description('<description>') on path '<path>' ]",
    description=f'Create a new workspace with an optional description and path.{info_workspaces}'
))

# Remove workspaces
statements.append(Forward(remove + workspace('workspace') + Word(alphas, alphanums + "_")('Workspace_Name'))('remove_workspace_statement'))
grammar_help.append(help_dict_create(
    name="remove workspace",
    category='Workspaces',
    command='remove workspace <workspace_name> ',
    description=f'Remove a workspace from your registry. Note that this doesn\'t remove the workspace\'s directory.{info_workspaces}'
))

# list workspaces
statements.append(Forward(lister + workspaces('workspaces'))('list_workspaces'))
grammar_help.append(help_dict_create(
    name="list workspaces",
    category='Workspaces',
    command="list workspaces",
    description=f'Lists all your workspaces.{info_workspaces}'
))

# endregion

##########################################################################
# region - Toolkits
# Note toolkits is the Caseless key word now .. simply changed in metadata
##########################################################################

info_toolkits_see_all = '\n<soft>To see all available toolkits, run <cmd>list all toolkits</cmd>.</soft>'
info_toolkits = '\n<soft>To learn more about toolkits, run <cmd>toolkit ?</cmd>.</soft>'


# Available commands per toolkit.
for tk in _all_toolkits:
    statements.append(Forward(CaselessKeyword(tk))(tk))
    grammar_help.append(help_dict_create(
        name=f'{tk} splash',
        category='Toolkits',
        command=tk.lower(),

        # This description is never read. Inside main.py -> do_help()
        # there is a clause that intercepts this command and displays
        # the available commands for the toolkit instead.
        description=f'Display the splash screen for the {tk} toolkit.'
    ))

# List toolkits
statements.append(Forward(lister + toolkits('toolkits'))('list_toolkits'))
grammar_help.append(help_dict_create(
    name="list toolkits",
    category='Toolkits',
    command="list toolkits",
    description=f'List all installed toolkits. To see all available toolkits, run <cmd>list all toolkits</cmd>.{info_toolkits_see_all}{info_toolkits}'
))

# List all toolkits
statements.append(Forward(lister + CaselessKeyword('all') + toolkits('toolkits'))('list_all_toolkits'))
grammar_help.append(help_dict_create(
    name="list all toolkits",
    category='Toolkits',
    command="list all toolkits",
    description=f'List all available toolkits.{info_toolkits}'
))

# Install toolkit
# Note: Update toolkit yet to be implemented (currently de-register, then add toolkit with new source)
statements.append(Forward(add + toolkit('workspace') + Word(alphas, alphanums + "_")('toolkit_name'))('add_toolkit'))
grammar_help.append(help_dict_create(
    name="add toolkit",
    category='Toolkits',
    command="add toolkit <toolkit_name>",
    description=f'Install a toolkit.{info_toolkits_see_all}{info_toolkits}'
))

# Remove toolkit
statements.append(Forward(remove + toolkit + Word(alphas, alphanums + "_")('toolkit_name'))('remove_toolkit'))
grammar_help.append(help_dict_create(
    name="remove toolkit",
    category='Toolkits',
    command="remove toolkit <toolkit_name>",
    description=('Remove a toolkit from the registry.\n'
                 'Note: This doesn\'t delete the toolkit code. If the toolkit is added again, a backup of the previous install is created in the toolkit directory at <yellow>~/.openad/toolkits</yellow>.'
                 f'{info_toolkits}')
    # Correct description but we have to update the functionality first.
    # description="Remove a toolkit from the registry. This affects all workspaces. A backup of the toolkit directory is stored in <yellow>~/.openad/toolkits_archive</yellow>.{info_toolkits}"
))

# Set a toolkit as the current context
statements.append(Forward(s_et + context('context') + Word(alphas, alphanums + "_")('toolkit_name') + Optional(CaselessKeyword('reset'))('reset'))('set_context'))
grammar_help.append(help_dict_create(
    name="set context",
    category='Toolkits',
    command="set context <toolkit_name> [ reset ]",
    description=f'Set your context to the chosen toolkit. By setting the context, the selected toolkit functions become available to you. The optional parameter \'reset\' can be used to reset your login information.{info_toolkits}'
))

# Get the current context
statements.append(Forward(get + context('context'))('get_context'))
grammar_help.append(help_dict_create(
    name="get context",
    category='Toolkits',
    command="get context",
    description='Display the currently selected toolkit.'
))

# Unset a toolkit as the current context
statements.append(Forward(unset + context('context'))('unset_context'))
grammar_help.append(help_dict_create(
    name="unset context",
    category='Toolkits',
    command="unset context",
    description=f'Exit your toolkit context. You will no longer have access to toolkit-specific functions.{info_toolkits}'
))

# endregion

##########################################################################
# region - Runs
##########################################################################

info_runs = '\n<soft>To learn more about runs, run <cmd>run ?</cmd>.</soft>'

# Create run
statements.append(Forward(create + run('run'))('create_run'))
grammar_help.append(help_dict_create(
    name="create run",
    category='Runs',
    command="create run",
    description=f'Start recording a run.{info_runs}'
))

# Save run
statements.append(Forward(save + run('run') + a_s + Word(alphas, alphanums + "_")('run_name'))('save_run'))
grammar_help.append(help_dict_create(
    name="save run",
    category='Runs',
    command="save run as <run_name>",
    description=f'Stop recording a run and save it.{info_runs}'
))

# Execute run
statements.append(Forward(run('run') + Word(alphas, alphanums + "_")('run_name'))('exec_run'))
grammar_help.append(help_dict_create(
    name="run",
    category='Runs',
    command="run <run_name>",
    description=f'Execute a previously recorded run. This will execute every command and continue regardless of any failures.{info_runs}'
))

# List runs
statements.append(Forward(lister + runs('runs'))('list_runs'))
grammar_help.append(help_dict_create(
    name="list runs",
    category='Runs',
    command="list runs",
    description=f'List all runs saved in the current workspace.{info_runs}'
))

# Display run
statements.append(Forward(d_isplay + run('run') + Word(alphas, alphanums + "_")('run_name'))('display_run'))
grammar_help.append(help_dict_create(
    name="display run",
    category='Runs',
    command="display run <run_name>",
    description=f'Display the commands stored in a certain run.{info_runs}'
))

# endregion

##########################################################################
# region - Utility
##########################################################################

# Display data
statements.append(Forward(d_isplay + data('data') + desc('file_path'))('display_data'))
grammar_help.append(help_dict_create(
    name="display data",
    category='Utility',
    command="display data '<csv_filename>'",
    description='Display data from a csv file.'
))

# --> result save --> Save data as csv
statements.append(Forward(result + save + Optional(a_s + desc('file_path')))('display_data__save'))
grammar_help.append(help_dict_create(
    name="save",
    category='Utility',
    command="result save [as '<csv_filename>']",
    description='Save table data to csv file.',
    parent="display data"
))

# --> result open --> Explore data in browser
statements.append(Forward(result + CaselessKeyword('open'))('display_data__open'))
grammar_help.append(help_dict_create(
    name="open",
    category='Utility',
    command="result open",
    description='Explore table data in the browser.',
    parent="display data"
))

# --> result edit --> Edit data in browser
statements.append(Forward(result + CaselessKeyword('edit'))('display_data__edit'))
grammar_help.append(help_dict_create(
    name="edit",
    category='Utility',
    command="result edit",
    description='Edit table data in the browser.',
    parent="display data"
))

# Edit config file (CLI-only)
if not is_notebook_mode():
    statements.append(Forward(CaselessKeyword('edit') + CaselessKeyword('config') + desc('json_file') + Optional(CaselessKeyword('schema') + desc('schema')))('edit_config'))
    grammar_help.append(help_dict_create(
        name="edit config",
        category='Utility',
        command="edit config '<json_config_file>' [ schema '<schema_file>']",
        description='Edit any JSON file in your workspace directly from the CLI. If a schema is specified, it will be used for validation and documentation.'
    ))

# Show molecules.
# Note: we don't allow dashes in dataframe names because it's a substraction operator and causes issues in Jupyter.
statements.append(Forward(
    show('show') + molecules + using + CaselessKeyword('dataframe') + Word(alphas, alphanums + "_")('in_dataframe')  # From dataframe
    + Optional(a_s + CaselessKeyword('molsobject')('object'))  # Return as molsobject
    + Optional(save + a_s + desc('results_file'))  # Save as csv/sdf
)('show_molecules_df'))
statements.append(Forward(
    show('show') + molecules + using + file + desc('moles_file')  # From mols file
    + Optional(a_s + CaselessKeyword('molsobject')('object'))  # Return as molsobject
    + Optional(save + a_s + desc('results_file'))  # Save as csv/sdf
)('show_molecules'))
grammar_help.append(help_dict_create(
    name="show molecules",
    category='Utility',
    command="show molecules using ( file '<mols_file>' | dataframe <dataframe> )\n    [ save as '<sdf_or_csv_file>' | as molsobject ]",
    description=f"""Launch the molecule viewer { 'in your browser ' if is_notebook_mode() else '' }to examine and select molecules from a SMILES sdf/csv dataset.

Examples:

    <cmd>show molecules using file 'base_molecules.sdf' as molsobject</cmd>
    <cmd>show molecules using dataframe my_dataframe save as 'selection.sdf'</cmd>
"""
))

# endregion

##########################################################################
# region - LLMs
##########################################################################

# Tell me chatbot
statements.append(Forward(CaselessKeyword('tell') + CaselessKeyword('me') + ZeroOrMore(Word(alphas, alphanums + "_" + '?' + "." + " " + "," + "'"))('Chat_String'))('how_do_i'))
grammar_help.append(help_dict_create(
    name="tell me",
    category='LLM',
    command="tell me <how to do xyz>",
    description='Ask your AI assistant how to do anything in OpenAD.'
))

statements.append(Forward(CaselessKeyword('set') + CaselessKeyword('llm') + ZeroOrMore(Word(alphas, alphanums + "_" + '?' + "." + " " + "," + "'"))('llm_name'))('set_llm'))
grammar_help.append(help_dict_create(
    name="set llm",
    category='LLM',
    command="set llm  <language_model_name>",
    description='Set the target language model name for the "tell me" command.'
))

statements.append(Forward(CaselessKeyword('clear') + CaselessKeyword('llm') + CaselessKeyword('auth'))('clear_llm_auth'))
grammar_help.append(help_dict_create(
    name="clear llm auth",
    category='LLM',
    command="clear llm auth",
    description='Clear the language model\'s authentication file.'
))


# endregion

##########################################################################
# region - File System
##########################################################################

# List files
statements.append(Forward(lister + CaselessKeyword('files'))('list_files'))
grammar_help.append(help_dict_create(
    name="list files",
    category='File System',
    command="list files",
    description='List all files in your current workspace.'
))

# Import file
statements.append(Forward(CaselessKeyword('import') + CaselessKeyword('from') + desc('source') + CaselessKeyword('to') + desc('destination'))('import_file'))
grammar_help.append(help_dict_create(
    name="import",
    category='File System',
    command="import from '<external_source_file>' to '<workspace_file>'",
    description='Import a file from outside OpenAD into your current workspace.'
))

# Export file
statements.append(Forward(CaselessKeyword('export') + CaselessKeyword('from') + desc('source') + CaselessKeyword('to') + desc('destination'))('export_file'))
grammar_help.append(help_dict_create(
    name="export",
    category='File System',
    command="export from '<workspace_file>' to '<external_file>'",
    description='Export a file from your current workspace to anywhere on your hard drive.'
))

# Copy file
statements.append(Forward(CaselessKeyword('copy') + CaselessKeyword('file') + desc('source') + CaselessKeyword('to') + desc('destination'))('copy_file'))
grammar_help.append(help_dict_create(
    name="copy",
    category='File System',
    command="copy file '<workspace_file>' to '<other_workspace_name>'",
    description='Export a file from your current workspace to another workspace.'
))

# Remove file
statements.append(Forward(CaselessKeyword('remove') + desc('file'))('remove_file'))
grammar_help.append(help_dict_create(
    name="remove",
    category='File System',
    command="remove '<filename>'",
    description='Remove a file from your current workspace.'
))

# endregion

##########################################################################
# region - Help
##########################################################################

# Display intro.
statements.append(Forward(CaselessKeyword('intro'))('intro'))
grammar_help.append(help_dict_create(
    name="intro",
    category='Help',
    command="intro",
    description='Display an introduction to the OpenAD CLI.'
))
# Open documentation webpage.
statements.append(Forward(CaselessKeyword('docs'))('docs'))
grammar_help.append(help_dict_create(
    name="docs",
    category='Help',
    command="docs",
    description='Open the documentation webpage.'
))
# List available commands
# Note - this is controlled directly from do_help.
grammar_help.append(help_dict_create(
    name="help",
    category='Help',
    command="?",
    description='List all available commands.'
))

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
grammar_help.append(help_dict_create(
    name="command help",
    category='Help',
    command="<soft>...</soft> ?",
    description='Display what a command does, or list all commands that contain this string.'
))

# endregion

##########################################################################
# region - Development
# The commands in this section are not intended for general use,
# and are not documented by the help system.
##########################################################################

# Launches the demo flask app.
statements.append(Forward(CaselessKeyword('flask') + CaselessKeyword('example'))('flask_example'))

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

    # if cmd_pointer.toolkit_current is None: # move down further to adress unset issue
    #    return
    global statements_zom

    cmd_pointer.current_statements_def = Forward()
    cmd_pointer.current_statements = orig_statements.copy()

    for i in orig_statements:
        cmd_pointer.current_statement_defs |= i
    if cmd_pointer.toolkit_current is not None:
        for i in cmd_pointer.toolkit_current.methods_grammar:
            cmd_pointer.current_statement_defs |= i
            cmd_pointer.current_statements.append(i)
    statements_zom = ZeroOrMore(statements_def)


def or_builder(options: list) -> str:
    expression = '('
    the_or = ''
    for i in options:
        expression = expression + the_or + i
        the_or = '|'

    return expression + ')'


def from_builder(options: list) -> str:
    clause = "CaselessKeyword('from')+"
    option_exp = []
    if not is_notebook_mode() and 'dataframe' in options:

        options.remove('dataframe')

    for i in options:
        if i == 'file':
            option_exp.append('(Suppress(CaselessKeyword("file"))+desc("from_file"))')
        elif i == 'dataframe':
            option_exp.append('(Suppress(CaselessKeyword("dataframe"))+' + 'Word(alphas, alphanums + "_")("from_dataframe"))')
        elif i == 'list':
            option_exp.append('(Suppress(CaselessKeyword("list"))+ Group(Suppress("[")+delimitedList(desc)("from_list")+Suppress("]")))')
    if len(option_exp) > 0:
        clause = clause + or_builder(option_exp)

    # elif len(option_exp) ==1:
    #    clause=clause+option_exp[0]
    else:
        raise ("invalid 'From' clause structure")

    return clause


def statement_builder(toolkit_pointer, statement):
    #####################################################################
    #
    # This section deals with Method call like statements, they will always be prefixed with 'exec'.
    ########################################################################################################
    if statement["exec_type"] == 'method':
        expression = "Suppress(e_xec)+" + 'CaselessKeyword ("' + statement['command'] + '")'
        expression = expression + '+ Suppress("(") +' + optional_parameter_list(statement, 'fixed_parameters') + ' +Suppress(")"))'

    elif statement["exec_type"] == 'standard_statement':
        key_leading_words = statement['command'].split()
        expression = ''
        plus = ''

        for i in key_leading_words:
            expression = expression + plus + 'CaselessKeyword ("' + i + '")'
            plus = '+'

        if 'SINGLE_PARM' in statement:
            if len(statement['SINGLE_PARM']) > 0:
                expression = expression + '+' + actual_parameter_list(statement, 'SINGLE_PARM')

        if 'from' in statement:
            if len(statement['from']) > 0:
                expression = expression + '+' + (from_builder(statement['from'])) + "('from_source')"

        if statement['USING'] is not None and len(statement['USING']) > 0:
            expression = expression + '+ Optional( (CaselessKeyword ("USING")+ Suppress("(") +' + optional_parameter_list(statement, 'USING') + '+Suppress(")") )("USING"))'

        if "RETURN_AS_DATA" in statement:
            expression = expression + '+ Optional(Suppress(CaselessKeyword ("return")) + Suppress(CaselessKeyword ("as"))+Suppress(CaselessKeyword ("data")))("return_as_data") '

        if "SAVE_AS" in statement:
            expression = expression + '+ Optional(Suppress(CaselessKeyword ("save")) + Suppress(CaselessKeyword ("as"))+desc("results_file") )("save_as") '
        if 'use_saved' in statement and len(statement['use_saved']) > 0:
            if statement['use_saved'] == 'True':
                expression = expression + '+' + "Optional(CaselessKeyword('use_saved'))('use_saved')"
        if 'async' in statement and len(statement['async']) > 0:
            if statement['async'] == 'both':
                expression = expression + '+' + "Optional(CaselessKeyword('async'))('do_async')"
            if statement['async'] == 'only':
                expression = expression + '+' + "CaselessKeyword('async')"

        expression = expression + ')'

    #####################################################################
    #
    # This section deals with search type statements.
    ########################################################################################################
    elif statement["exec_type"] == 'search_statement':
        expression = ""
        a_char = ' '
        for i in statement['command'].split():

            expression = expression + a_char + ' CaselessKeyword ("' + i + '")'
            if a_char == ' ':
                a_char = ' + '
        expression = expression + " +desc('val')('" + statement['subject'] + "') "
        ####################################################################################################
        #
        # "USING" clase specifies any indexes, attributes or other services that optionally can be used
        # currently if a field is mandatory this will need to be picked up in the toolkit program
        # e.g. Using index_type="arxiv" contstraint=chemicals maxhops=1
        ####################################################################################################
        if statement['USING'] is not None:
            expression = expression + '+ (CaselessKeyword ("USING")+ Suppress("(") +' + \
                optional_parameter_list(
                    statement, 'USING') + '+Suppress(")") )("USING")'
        ######################################################################################################
        #
        # Show (<data_types>) is to specify what data to retrieve.
        # This is optional clause, any mandatory values need to be managed by the toolkit.
        # ####################################################################################################

        if "SHOW" in statement:
            if statement['SHOW'] is not None:
                options = []
                for i in statement['SHOW']:
                    options.append(i)

                expression = expression + '+ Optional(CaselessKeyword ("SHOW") +Suppress("(")+' + \
                    oneormore_str(options) + ' +Suppress(")"))("show_data") '
        ######################################################################################################
        #
        # Estimate Only is for statements that may take a long time and there is an estimate capability.
        # ####################################################################################################
        if "ESTIMATE_ONLY" in statement:
            expression = expression + \
                '+ Optional(Suppress(CaselessKeyword ("ESTIMATE")) +Suppress(CaselessKeyword ("ONLY")) )("estimate_only") '
        if "RETURN_AS_DATA" in statement:
            expression = expression + '+ Optional(Suppress(CaselessKeyword ("return")) + Suppress(CaselessKeyword ("as"))+Suppress(CaselessKeyword ("data")))("return_as_data") '

        if "SAVE_AS" in statement:
            expression = expression + \
                '+ Optional(Suppress(CaselessKeyword ("save")) + Suppress(CaselessKeyword ("as"))+desc("results_file") )("save_as") '

        expression = expression + ')'

    try:

        toolkit_pointer.methods_grammar.append(eval(" Forward( " + expression + ' ("toolkit_exec_' + statement['command'] + '")'))

        toolkit_pointer.methods.append(statement['command'])

        toolkit_pointer.methods_execute.append(statement['method'])

        toolkit_pointer.methods_library.append(statement['library'])

        toolkit_pointer.methods_dict.append(statement)

        # Clause Amendment
        clause_Amendment = ''
        if 'USING' in statement:
            if statement['USING'] is not None:
                clause_Amendment = clause_Amendment + '\n\n NOTE: The Using Clause Requires all the Parameters added to the Using Clause be in the defined order as per in the above help documentation'

        if clause_Amendment != '':

            statement['help']['description'] = statement['help']['description'] + clause_Amendment

        toolkit_pointer.methods_help.append(statement['help'])

    except BaseException as err:
        fwd_expr = 'Forward( ' + expression + ' ("toolkit_exec_' + statement['command'] + '")'
        output_error(msg('err_add_command', statement['command'], fwd_expr, err, split=True))

    return True


def oneormore_str(options: list):
    expression = ' OneOrMore('
    divider = ' '
    ii = 0
    for i in options:
        expression = expression + divider + \
            ' CaselessKeyword ("' + \
            i.strip() + '") '

        divider = "|"

    expression = expression + ")"
    return expression


def optional_parameter_list(statement: dict, clause: str):
    ii = 0
    expression = " "
    for i in statement[clause]:
        if ii == 0:
            expression = expression + " "
            if statement[clause][i] == 'str':
                expression = expression + \
                    " ZeroOrMore(Group( CaselessKeyword ('" + i + \
                    "') +Suppress('=')+key_val_expr('val'))('" + \
                    i + "'))" + " "
            elif statement[clause][i] == 'desc':
                expression = expression + \
                    " ZeroOrMore(Group( CaselessKeyword ('" + i + \
                    "') +Suppress('=')+desc('val'))('" + \
                    i + "'))" + " "
            else:
                expression = expression + \
                    " ZeroOrMore(Group( CaselessKeyword ('" + i + \
                    "') +Suppress('=')+key_val_expr_num('val'))('" + \
                    i + "'))" + " "
        else:

            if statement[clause][i] == 'str':
                expression = expression + \
                    " & ZeroOrMore(Group( CaselessKeyword ('" + i + \
                    "') +Suppress('=')+key_val_expr('val'))('" + \
                    i + "'))" + " "
            elif statement[clause][i] == 'desc':
                expression = expression + \
                    " & ZeroOrMore(Group( CaselessKeyword ('" + i + \
                    "') +Suppress('=')+desc('val'))('" + \
                    i + "'))" + " "
            else:
                expression = expression + \
                    " & ZeroOrMore(Group( CaselessKeyword ('" + i + \
                    "') +Suppress('=')+key_val_expr_num('val'))('" + \
                    i + "'))" + " "
        ii = 1

    return expression


def actual_parameter_list(statement: dict, clause: str):
    ii = 0
    expression = " "
    for i in statement[clause]:
        if ii == 0:
            expression = expression + " "
            if statement[clause][i] == 'str':
                expression = expression + \
                    " key_val_expr('" + i + "') "
            elif statement[clause][i] == 'desc':
                expression = expression + \
                    "desc('" + i + "') "
            else:
                expression = expression + \
                    " key_val_expr_num('" + i + "') "
        else:

            if statement[clause][i] == 'str':
                expression = expression + \
                    ", key_val_expr('" + i + "') "
            elif statement[clause][i] == 'desc':
                expression = expression + \
                    ",desc('" + i + "') "
            else:
                expression = expression + \
                    ", key_val_expr_num('" + i + "') "
        ii = 1
    return expression


def output_train_statements(cmd_pointer):
    import os
    import glob
    training_statements = []
    toolkit_lists = []
    toolkit_name = []
    toolkit_description = []
    i = 0
    cmd_pointer.home_dir

    try:
        if not os.path.exists(os.path.expanduser(os.path.expanduser(cmd_pointer.home_dir + "/prompt_train/"))):
            os.mkdir(os.path.expanduser(os.path.expanduser(cmd_pointer.home_dir + "/prompt_train/")))
    except BaseException:
        raise Exception(f"unable to create prompt soc directory {str(os.path.expanduser(cmd_pointer.home_dir+'/prompt_train/'))}")

    for file in glob.glob(os.path.expanduser(str(os.path.expanduser(cmd_pointer.home_dir + '/prompt_train/')) + "/*")):
        os.remove(file)

    # while i < len(orig_statements): @Phil this causes issues when statements and help are out of sync, eg. for display molecules.
    while i < len(grammar_help):
        training_statements.append({'command_group': 'base', 'command_name': grammar_help[i]['name'], 'command_syntax': grammar_help[i]['command'], 'command_help': grammar_help[i]['description']})
        i += 1

    file = open(os.path.expanduser(cmd_pointer.home_dir + '/prompt_train/base_commands.cdoc'), 'w', newline='\n')
    file.write("""openad client information

        For information in providing answers to how to or Help questions from users :



        The below describes cdccl clients domain specific language (DSL) for managing science activities using the DSL

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
            Display: display a file or result set
            Show:  Show a data set using a utility that enables you to manipulate or diagramatically view it.
            Backup: backup a plugin or workspace
            Add: add a function or plugin
            Remove: delete an object
            Save: Save a run or file of some kind
            Load: load a file from project directory to Target system
            pyparsing_statement: a statement defined using pyparsing for the domain specific language
            help_text: description of the Domain Specific language statement defined in a pyparsing_statement
            toolkit: these are contextual plugins that are available one at a time for providing specific functionality to the user. Valid toolkits are DS4SD (deepSearch), GT4SD(generative AI toolkit), RXN (retro synthesis), ST4SD(simulation toolkit)
            History: History of DSL commands for a given Workspace
            run: list of sequential commands saved by the user')
            ?: will display help and if positioned prior to a command will display help options for that command \\@ \n\n""")

    file.write('The following dictionaries are the Domain Specific Language (DSL) commands available in the base openad client.\n  The users Domain Specific Language commands will be interpreted from these command definitions which the application trainslates into the Domain specific Language (DSL) .. for reference \n command_name : is the name of the command\n command_group : is the toolkit group the command belongs to, command_syntax : is the Syntax description for the command \n command_help : is the associated help instructions and examples on how to use the command.  \\@\n')
    for i in training_statements:
        file.write(str(i) + '\\@\n')
    file.close()

    for i in cmd_pointer.settings['toolkits']:
        token, a_toolkit = load_toolkit(i)
        training_statements = []
        file = open(os.path.expanduser(cmd_pointer.home_dir + '/prompt_train/toolkit_' + i + '.cdoc'), 'w', newline='\n')
        file.write('This file represents the description of the ' + a_toolkit.toolkit_name + ' \n\n')

        file.write('this is the description of the overall capability of the ' + a_toolkit.toolkit_name + ' \n \n')
        file.write('The following are the  command definitions for the ' + a_toolkit.toolkit_name + ' toolkit with their  command definition and description \n ' +
                   str(a_toolkit.toolkit_description) + ' \n When displaying an answer always interpret the command specified \\@ \n')

        x = 0
        while x < len(a_toolkit.methods):
            training_statements.append({'toolkit group': a_toolkit.toolkit_name, 'command_name': a_toolkit.methods_help[x]['name'], 'command': a_toolkit.methods_help[x]['command'], 'command_help': a_toolkit.methods_help[x]['description']})
            x += 1
        for i in training_statements:
            file.write(str(i) + '\n\\@',)
        file.close()
    return


# Need to fix later
#!/usr/local/opt/python@3.9/bin/python3.9
import os
import imp
import glob
import json
import logging
from ad4e_opentoolkit.helpers.output import msg, output_error

# Globals
from ad4e_opentoolkit.app.global_var_lib import _meta_dir_toolkits
from ad4e_opentoolkit.app.global_var_lib import _meta_workspaces


class Toolkit:
    def __init__(self, name) -> None:
        self.toolkit_name = name
        self.toolkit_description = None
        self.methods = []
        self.methods_grammar = []
        self.methods_execute = []
        self.methods_help = []
        self.methods_dict = []
        self.methods_library = []


# Load all toolkit statments.
def load_toolkit(toolkit_name):
    the_toolkit = Toolkit(toolkit_name)

    for i in glob.glob(_meta_dir_toolkits + "/" + toolkit_name + "/**/func_*.json", recursive=True):
        func_file = open(i, "r")
        x = json.load(func_file)
        statement_builder(the_toolkit, x)
    try:
        with open(_meta_dir_toolkits + "/" + toolkit_name + "/description.txt", 'r')as file:
            the_toolkit.toolkit_description = file.read()
            file.close()
    except BaseException:
        # print('no description')
        the_toolkit.toolkit_description = None
    return True, the_toolkit
