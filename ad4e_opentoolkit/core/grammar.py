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
import ad4e_opentoolkit.core.help as adccl_help

# Helpers
from ad4e_opentoolkit.helpers.general import is_notebook_mode
from ad4e_opentoolkit.helpers.output import output_error, msg

# Global variables
from ad4e_opentoolkit.app.global_var_lib import _all_toolkits

# @Phil can be deleted ?
# from pyparsing import *
# from operator import ior


get, lister, description, using, create, s_et, unset, workspace, workspaces, context, jobs, e_xec, a_s, optimize, w_ith, toolkits, toolkit, GPU, experiment, add, run, save, runs, show, molecules, file, d_isplay, history, data, remove = \
    map(CaselessKeyword, "get list description using create set unset workspace workspaces context jobs exec as optimize with toolkits toolkit gpu experiment add run save runs show molecules file display history data remove".split())
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

statements = []  # Statement Definitions stored here
grammar_help = []  # Help text stored here


##########################################################################
# region - General Vocabulary
##########################################################################

# General splash screen
statements.append(Forward(CaselessKeyword('adccl'))('welcome'))
grammar_help.append(adccl_help.help_dict_create(
    name="welcome",
    category='General',
    command="adccl",
    description='Display the ADCCL splash screen'
))

# Get status
statements.append(Forward(get + CaselessKeyword('status'))('get_status'))
grammar_help.append(adccl_help.help_dict_create(
    name="get status",
    category='General',
    command="get status",
    description='Return the current workspace and toolkit for the Utility'
))

# Display history
statements.append(Forward(d_isplay + history('history'))('display_history'))
grammar_help.append(adccl_help.help_dict_create(
    name="Display History",
    category='General',
    command="display history ",
    description='This command displays the last 30 commands run in the Workspace'
))

# Display data
statements.append(Forward(d_isplay + data('data') + desc('file_path'))('display_data'))
grammar_help.append(adccl_help.help_dict_create(
    name="display data",
    category='General',
    command="display data '<csv_filename>'",
    description='Display the data from a csv file.'
))

# --> save --> Save data as csv
statements.append(Forward(save + Optional(CaselessKeyword('as') + desc('file_path')))('display_data__save'))
grammar_help.append(adccl_help.help_dict_create(
    name="save",
    category='General',
    command="save [as '<csv_filename>']",
    description='Store table data into a csv file.'
))

# --> open --> Explore data in browser
statements.append(Forward(CaselessKeyword('open'))('display_data__open'))
grammar_help.append(adccl_help.help_dict_create(
    name="save",
    category='General',
    command="save [as '<csv_filename>']",
    description='Store table data into a csv file.'
))

# Clear sessions
statements.append(Forward(CaselessKeyword('clear') + CaselessKeyword('sessions'))('clear_sessions'))
grammar_help.append(adccl_help.help_dict_create(
    name="Clear Sessions",
    category='General',
    command="clear sessions",
    description='Clears any other sessions that may be running.'
))

# How do I chat
statements.append(Forward(CaselessKeyword('tell') + CaselessKeyword('me') + ZeroOrMore(Word(alphas, alphanums + "_" + '?' + "." + " " + "," + "'"))('Chat_String'))('HOW_DO_I'))
grammar_help.append(adccl_help.help_dict_create(
    name="tell me",
    category='LLM',
    command="tell me  <chat string>",
    description='This command asks chatgpt how to do set of tasks in the adccl client '
))

statements.append(Forward(CaselessKeyword('set') + CaselessKeyword('llm') + ZeroOrMore(Word(alphas, alphanums + "_" + '?' + "." + " " + "," + "'"))('llm_name'))('set_llm'))
grammar_help.append(adccl_help.help_dict_create(
    name="set llm",
    category='LLM',
    command="set llm  <name of language model>",
    description='This command sets the target language model type for the "tell me" command '
))

statements.append(Forward(CaselessKeyword('clear') + CaselessKeyword('llm') + CaselessKeyword('auth'))('clear_llm_auth'))
grammar_help.append(adccl_help.help_dict_create(
    name="clear llm auth",
    category='LLM',
    command="clear llm auth",
    description='This command clears the llm authority file. '
))


# endregion

##########################################################################
# region - Workspace Vocabulary
##########################################################################

# Wildcard - this lets us catch any unrecognized commands.
# statements.append(Forward(Word(string.printable))('_'))

# Set workspaces
statements.append(Forward(s_et + workspace('workspace') + Word(alphas, alphanums + "_")('Workspace_Name'))('set_workspace_statement'))
grammar_help.append(adccl_help.help_dict_create(
    name="set workspace",
    category='Workspaces',
    command="set workspace <workspace_name>",
    description='Change the current workspace.'
))

# Create workspaces
statements.append(Forward(create + workspace('workspace') + Word(alphas, alphanums + "_")('Workspace_Name') +
                          Optional(description('description') + Suppress("(") + desc("proj_desc") + Suppress(")"))('description_option') +
                          Optional(CaselessKeyword('on') + CaselessKeyword('path') + desc("w_path"))('workspace_path'))('create_workspace_statement'))
grammar_help.append(adccl_help.help_dict_create(
    name="create workspace",
    category='Workspaces',
    command="create workspace <workspace_name> [ description('<workspace_description>') on path '<path>' ]",
    description=('Create a new workspace with an optional description. '
                 'To learn more about what a workspace is, run <cmd>workspace ?</cmd>')
))

# Remove workspaces
statements.append(Forward(remove + workspace('workspace') + Word(alphas, alphanums + "_")('Workspace_Name'))('remove_workspace_statement'))
grammar_help.append(adccl_help.help_dict_create(
    name="remove workspace",
    category='Workspaces',
    command='remove workspace <workspace_name> ',
    description=('Remove a workspace from your registry. This doesn\'t remove the workspace\'s directory '
                 'from the adccl home directory ($HOME/.addcl/workspaces). '
                 'To learn more about what a workspace is, run <cmd>workspace ?</cmd>')
))

# list workspaces
statements.append(Forward(lister + workspaces('workspaces'))('list_workspaces'))
grammar_help.append(adccl_help.help_dict_create(
    name="list workspaces",
    category='Workspaces',
    command="list workspaces",
    description='This command lists all workspaces '
))

# Get workspace
statements.append(Forward(get + workspace('workspace') + Word(alphas, alphanums + "_")('Workspace_Name'))('get_workspace'))
grammar_help.append(adccl_help.help_dict_create(
    name="get workspace",
    category='Workspaces',
    command="get workspace <workspace_name>",
    description='This command retrieves details about the workspace '
))

# endregion

##########################################################################
# region - Toolkit Vocabulary
# Note toolkits is the Caseless key word now .. simply changed in metadata
##########################################################################

# Splash screens for toolkits
for tk in _all_toolkits:
    statements.append(Forward(CaselessKeyword(tk))(tk))
    grammar_help.append(adccl_help.help_dict_create(
        name=f'{tk} splash',
        category='Toolkits',
        command=tk.lower(),

        # We display the full toolkit help when you type `<toolkit_name> ?`
        description=f'Display the splash screen for the {tk} toolkit.'
    ))

# Remove toolkit
statements.append(Forward(remove + toolkit + Word(alphas, alphanums + "_")('toolkit_name'))('remove_toolkit'))
grammar_help.append(adccl_help.help_dict_create(
    name="remove toolkit",
    category='Toolkits',
    command="remove toolkit <toolkit_name>",
    # description=('De-register a toolkit from the registry. This doesn\'t delete the toolkit code. '
    #              'If the toolkit is added again, a backup of the directories containing the code is placed in the main toolkit directory.')

    # REVIEW
    description="Remove a toolkit from the registry. This affects all workspaces. A backup of the toolkit directory is stored in '~/.adccl/toolkits_archive'."

))

# Install toolkit
# Note: Update toolkit yet to be implemented (currently de-register, then add toolkit with new source)
statements.append(Forward(add + toolkit('workspace') + Word(alphas, alphanums + "_")('toolkit_name'))('add_toolkit'))
grammar_help.append(adccl_help.help_dict_create(
    name="add toolkit",
    category='Toolkits',
    command="add toolkit <toolkit_name>",
    description='Install a toolkit. To see all available toolkits, run <cmd>list all toolkits</cmd>.'
))

# List toolkits
statements.append(Forward(lister + toolkits('toolkits'))('list_toolkits'))
grammar_help.append(adccl_help.help_dict_create(
    name="list toolkits",
    category='Toolkits',
    command="list toolkits",
    description='Display a list of installed toolkits. To see all available toolkits, run <cmd>list all toolkits</cmd>.'
))

# List all toolkits
statements.append(Forward(lister + CaselessKeyword('all') + toolkits('toolkits'))('list_all_toolkits'))
grammar_help.append(adccl_help.help_dict_create(
    name="list all toolkits",
    category='Toolkits',
    command="list all toolkits",
    description='Display a list of all available toolkits.'
))

# Set a toolkit as the current context
statements.append(Forward(s_et + context('context') + Word(alphas, alphanums + "_")('toolkit_name') + Optional(CaselessKeyword('reset'))('reset'))('set_context'))
grammar_help.append(adccl_help.help_dict_create(
    name="set context",
    category='Toolkits',
    command="set context <toolkit_name> [reset]",
    description='Set your function context to the chosen toolkit. The context dictates functions available to you. The optional parameter reset can be used to reset login information'
))

# Unset a toolkit as the current context
statements.append(Forward(unset + context('context'))('unset_context'))
grammar_help.append(adccl_help.help_dict_create(
    name="unset context",
    category='Toolkits',
    command="unset context",
    description='Exit your toolkit context. You will no longer have access to toolkit-specific functions.'
))

# endregion

##########################################################################
# region - Run Vocabulary
##########################################################################

# Create run
statements.append(Forward(create + run('run'))('create_run'))
grammar_help.append(adccl_help.help_dict_create(
    name="create run",
    category='Runs',
    command="create run",
    description='This command flags a position for the beginning of a run recording.'
))

# Save run
statements.append(Forward(save + run('run') + a_s + Word(alphas, alphanums + "_")('run_name'))('save_run'))
grammar_help.append(adccl_help.help_dict_create(
    name="save run",
    category='Runs',
    command="save run as <run_name>",
    description='This command saves a run up until the last point in time a "create run" has been called'
))

# Execute run
statements.append(Forward(run('run') + Word(alphas, alphanums + "_")('run_name'))('exec_run'))
grammar_help.append(adccl_help.help_dict_create(
    name="run",
    category='Runs',
    command="run <run_name>",
    description='This command iterates through a previously recorded run. It will execute every command and continue whether there is a failure or not'
))

# List runs
statements.append(Forward(lister + runs('runs'))('list_runs'))
grammar_help.append(adccl_help.help_dict_create(
    name="list runs",
    category='Runs',
    command="list runs",
    description='This command lists all runs for the current workspace'
))

# Display run
statements.append(Forward(d_isplay + run('run') + Word(alphas, alphanums + "_")('run_name'))('display_run'))
grammar_help.append(adccl_help.help_dict_create(
    name="display run",
    category='Runs',
    command="display run <run_name>",
    description='This command displays the contents of a run files'
))

# endregion

##########################################################################
# region - File System Vocabulary
##########################################################################

# List files
statements.append(Forward(lister + CaselessKeyword('files'))('list_files'))
grammar_help.append(adccl_help.help_dict_create(
    name="list files",
    category='File System',
    command="list files",
    description='List all files in your current workspace.'
))

# Import file
statements.append(
    Forward(
        CaselessKeyword('import') +
        CaselessKeyword('from') +
        desc('source') +
        CaselessKeyword('to') +
        desc('destination'))('import_file'))
grammar_help.append(adccl_help.help_dict_create(
    name="import",
    category='File System',
    command="import from '<external_source_file>' to '<workspace_file>'",
    description='Import files from outside adccl into the current workspace'
))

# Export file
statements.append(
    Forward(
        CaselessKeyword('export') +
        CaselessKeyword('from') +
        desc('source') +
        CaselessKeyword('to') +
        desc('destination'))('export_file'))
grammar_help.append(adccl_help.help_dict_create(
    name="export",
    category='File System',
    command="export from '<workspace_file>' to '<external_file>'",
    description='exports a file from current workspace to an external file'
))

# Copy file
statements.append(
    Forward(
        CaselessKeyword('copy') +
        CaselessKeyword('file') +
        desc('source') +
        CaselessKeyword('to') +
        desc('destination'))('copy_file'))
grammar_help.append(adccl_help.help_dict_create(
    name="copy",
    category='File System',
    command="copy from '<workspace_file>' to '<other_workspace_name>'",
    description='exports a file from current workspace to another workspace'
))

# Remove file
statements.append(Forward(CaselessKeyword('remove') + desc('file'))('remove_file'))
grammar_help.append(adccl_help.help_dict_create(
    name="rm",
    category='File System',
    command="remove '<filename>'",
    description='remove files from workspaces'
))

# endregion


# Define The Concepts of Jobs
# statements.append( Forward(lister +jobs('job') )('list_jobs'))

# This is a space that will require further thought whether to run jobs
# asynchronously from the client and how to track and save their results


##########################################################################
# region - Global viewers
##########################################################################

if is_notebook_mode():
    # Jupyter

    # Show molecules using file.
    statements.append(Forward(show('show') + molecules + using + file + desc('moles_file') + Optional(CaselessKeyword('as') + CaselessKeyword('molsobject')('object')))('show_molecules'))
    grammar_help.append(adccl_help.help_dict_create(
        name="show molecules",
        category='GUI',
        command="show molecules using file <mols_file> [  as molsobject ]",
        description='This command shows 1 or more Molecules in selection grid using a dataset containing smiles molecule strings. It display a molecule selection grid inside Jypyter Notebooks or Jupyter Lab. \n example command:  "show molecules using dataframe mydata_frame as molsobject "  '
    ))

    # Show molecules using dataframe.
    statements.append(Forward(show('show') + molecules + using + CaselessKeyword('dataframe') + Word(alphas, alphanums + "_")('in_dataframe') +
                      Optional(CaselessKeyword('as') + CaselessKeyword('molsobject')('object')))('show_api_molecules'))
    grammar_help.append(adccl_help.help_dict_create(
        name="show molecules",
        category='GUI',
        command="show molecules using dataframe <dataframe> [  as molsobject ]",
        description='This command shows 1 or more Molecules in selection grid using a dataset containing smiles molecule strings. It display a molecule selection grid inside Jypyter Notebooks or Jupyter Lab. \n example command:  "show molecules using file \'base_molecules.sdf\' as molsobject "'
    ))
else:
    # CLI

    # Show molecules using file.
    statements.append(Forward(show('show') + molecules + using + file + desc('moles_file') + Optional(
        CaselessKeyword('Save') + CaselessKeyword('as') + desc('results_file')))('show_molecules'))
    grammar_help.append(adccl_help.help_dict_create(
        name="show molecules",
        category='GUI',
        command="show molecules using file '<mols_file>' [ save as '<sdf_or_csv_file>' ]",
        description='''This command shows 1 or more Molecules in selection grid using a dataset containing smiles molecule strings. It is displayed via browser when run from a command line and only accepts files of type .sdf or .csv.
        Input files may only be of .type sdf or .csv and you can only save files as .sdf or .csv types.
        Files referenced are from your local workspace, it does not support referencing files outside your workspace. If you wish to move a file into your workspace see the import command.
        Examples:
               - "show molecules using file 'my_molecules.sdf'"
               - "show molecules using file 'my_molecules.sdf' Save as 'selection.csv' "

        Alternatively when run in jupyter notebooks can accept a dataframe. See help when in Notebooks for details. '''
    ))

    # Edit config file.
    statements.append(Forward(CaselessKeyword('edit') + CaselessKeyword('config') + desc('json_file') + Optional(CaselessKeyword('template') + desc('template')))('edit_config'))
    grammar_help.append(adccl_help.help_dict_create(
        name="edit config",
        category='General',
        command="edit config '<json_config_file>' [template '<template_file>']",
        description='edits a config card in a workspace, if template specified will use template. All files assumed to be in the current workspace directory.'
    ))

# endregion

##########################################################################
# region - Help Vocabulary
##########################################################################

# Display intro.
statements.append(Forward(CaselessKeyword('intro'))('intro'))
grammar_help.append(adccl_help.help_dict_create(
    name="intro",
    category='Help',
    command="intro",
    description='Displays Introduction Help Text'
))
# Open documentation webpage.
statements.append(Forward(CaselessKeyword('docs'))('docs'))
grammar_help.append(adccl_help.help_dict_create(
    name="docs",
    category='Help',
    command="docs",
    description='Open the documentation webpage.'
))
# List available commands
# Note - this is controlled directly from do_help.
grammar_help.append(adccl_help.help_dict_create(
    name="help",
    category='Help',
    command="?",
    description='Lists all available commands.'
))

# Open advanced help.
# Note - this is controlled directly from do_help.
grammar_help.append(adccl_help.help_dict_create(
    name="advanced help",
    category='Help',
    command="??",
    description='Displays interactive help interface in the CLI.'
))

# Display command help
# Note - this is controlled directly from do_help.
grammar_help.append(adccl_help.help_dict_create(
    name="command help",
    category='Help',
    command="<soft>...</soft> ?",
    description='Explains what a command does, or lists all commands that either contain this word or start with this string.'
))

# endregion

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
    if cmd_pointer.toolkit_current is None:
        return
    global statements_zom

    cmd_pointer.current_statements_def = Forward()
    cmd_pointer.current_statements = orig_statements.copy()
    for i in orig_statements:
        cmd_pointer.current_statement_defs |= i
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
                    " Optional(Group( CaselessKeyword ('" + i + \
                    "') +Suppress('=')+key_val_expr('val'))('" + \
                    i + "'))" + " "
            elif statement[clause][i] == 'desc':
                expression = expression + \
                    " Optional(Group( CaselessKeyword ('" + i + \
                    "') +Suppress('=')+desc('val'))('" + \
                    i + "'))" + " "
            else:
                expression = expression + \
                    " Optional(Group( CaselessKeyword ('" + i + \
                    "') +Suppress('=')+key_val_expr_num('val'))('" + \
                    i + "'))" + " "
        else:

            if statement[clause][i] == 'str':
                expression = expression + \
                    "+ Optional(Group( CaselessKeyword ('" + i + \
                    "') +Suppress('=')+key_val_expr('val'))('" + \
                    i + "'))" + " "
            elif statement[clause][i] == 'desc':
                expression = expression + \
                    "+ Optional(Group( CaselessKeyword ('" + i + \
                    "') +Suppress('=')+desc('val'))('" + \
                    i + "'))" + " "
            else:
                expression = expression + \
                    "+ Optional(Group( CaselessKeyword ('" + i + \
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

    while i < len(orig_statements):
        training_statements.append({'command_group': 'base', 'command_name': grammar_help[i]['name'], 'command_syntax': grammar_help[i]['command'], 'command_help': grammar_help[i]['description']})
        i += 1
    file = open(os.path.expanduser(cmd_pointer.home_dir + '/prompt_train/base_commands.cdoc'), 'w', newline='\n')
    file.write("""adccl client information

        For information in providing answers to how to or Help questions from users :



        The below describes cdccl clients domain specific language (DSL) for managing science activities using the DSL

        Vocabulary:
            DSL: Domain Specific Language or DSL  that is implemented for the ADCCL client and all commands are formatted in
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

    file.write('The following dictionaries are the Domain Specific Language (DSL) commands available in the base adccl client.\n  The users Domain Specific Language commands will be interpreted from these command definitions which the application trainslates into the Domain specific Language (DSL) .. for reference \n command_name : is the name of the command\n command_group : is the toolkit group the command belongs to, command_syntax : is the Syntax description for the command \n command_help : is the associated help instructions and examples on how to use the command.  \\@\n')
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
from ad4e_opentoolkit.core.grammar import statement_builder
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
