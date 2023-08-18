#!/usr/local/opt/python@3.9/bin/python3.9
import os
import ad4e_opentoolkit.app.login_manager as login_manager
import readline

# Core
from ad4e_opentoolkit.core.lang_file_system import import_file, export_file, copy_file, remove_file, list_files
from ad4e_opentoolkit.core.lang_sessions_and_registry import clear_other_sessions, write_registry, registry_add_toolkit, registry_deregister_toolkit, initialise_registry,update_main_registry_env_var
from ad4e_opentoolkit.core.lang_workspaces import create_workspace, remove_workspace, list_workspaces, set_workspace, get_workspace
from ad4e_opentoolkit.core.lang_mols import display_mols
from ad4e_opentoolkit.core.grammar import create_statements
from ad4e_opentoolkit.core.lang_runs import display_run, execute_run, save_run, list_runs

# Toolkits
from ad4e_opentoolkit.toolkit.toolkit_main import load_toolkit
from ad4e_opentoolkit.toolkit.toolkit_main import execute_tookit
from ad4e_opentoolkit.toolkit.toolkit_main import load_toolkit_description
from  ad4e_opentoolkit.llm_assist.llm_interface import how_do_i,set_llm,clear_llm_auth

# Global variables
from ad4e_opentoolkit.app.global_var_lib import _meta_dir as _meta_dir
from ad4e_opentoolkit.app.global_var_lib import _meta_dir_toolkits as _meta_dir_toolkits
from ad4e_opentoolkit.app.global_var_lib import _meta_registry as _meta_registry
from ad4e_opentoolkit.app.global_var_lib import _meta_registry_session as _meta_registry_session
from ad4e_opentoolkit.app.global_var_lib import _meta_login_registry as _meta_login_registry
from ad4e_opentoolkit.app.global_var_lib import _meta_workspaces as _meta_workspaces
from ad4e_opentoolkit.app.global_var_lib import _all_toolkits as _all_toolkits

# Helpers
from ad4e_opentoolkit.helpers.output import msg, output_text, output_error, output_warning, output_success, output_table
from ad4e_opentoolkit.helpers.style_parser import tags_to_markdown
from ad4e_opentoolkit.helpers.general import refresh_prompt
from ad4e_opentoolkit.helpers.splash import splash
from ad4e_opentoolkit.helpers.output_content import adccl_intro


# This is called by the default run_cmd method, for executing current commands.
def lang_parse(cmd_pointer, parser):
    
    # Workspace commands
    if parser.getName() == 'create_workspace_statement':
        return create_workspace(cmd_pointer, parser)  # Addressed
    elif parser.getName() == 'remove_workspace_statement':
        return remove_workspace(cmd_pointer, parser)  # Addressed
    elif parser.getName() == 'set_workspace_statement':
        return set_workspace(cmd_pointer, parser)
    elif parser.getName() == 'list_workspaces':
        return list_workspaces(cmd_pointer, parser)
    elif parser.getName() == 'get_workspace':
        return get_workspace(cmd_pointer, parser)

    # Toolkit welcome screens
    elif parser.getName() in _all_toolkits:
        return output_text(splash(parser.getName(), cmd_pointer, raw=True), nowrap=True)

    # Toolkit commands
    elif parser.getName() == 'add_toolkit':
        
        return registry_add_toolkit(cmd_pointer, parser)
    elif parser.getName() == 'remove_toolkit':
        return registry_deregister_toolkit(cmd_pointer, parser)
    elif parser.getName() == 'list_toolkits':
        return list_toolkits(cmd_pointer, parser)
    elif parser.getName() == 'list_all_toolkits':
        return list_all_toolkits(cmd_pointer, parser)
    elif parser.getName() == 'set_context':
        return set_context(cmd_pointer, parser)
    elif parser.getName() == 'unset_context':
        return unset_context(cmd_pointer, parser)
    
    #Language Model How To
    elif parser.getName() == 'HOW_DO_I':
        result= how_do_i(cmd_pointer,parser)
        if result == False:
            
            return False
        update_main_registry_env_var(cmd_pointer,'refresh_help_ai',False)
        write_registry(cmd_pointer.settings,cmd_pointer)
        return result
    elif parser.getName() == 'set_llm':
        try:
            result= set_llm(cmd_pointer,parser)          
            cmd_pointer.llm_model=cmd_pointer.llm_models[cmd_pointer.llm_service] 
            update_main_registry_env_var(cmd_pointer,'llm_service',cmd_pointer.llm_service)
            cmd_pointer.refresh_vector=True
            cmd_pointer.refresh_train=True
            write_registry(cmd_pointer.settings,cmd_pointer,False)
            write_registry(cmd_pointer.settings,cmd_pointer,True)
        except Exception as e:
            print(e)
        return result
    elif parser.getName() == 'clear_llm_auth':
        result= clear_llm_auth(cmd_pointer,parser)
        return result

    # Run commands
    elif parser.getName() == 'create_run':
        # This simply works off the history file entry no procedure needed.
        return output_text(msg('create_run_started'), cmd_pointer, pad=1, nowrap=True)
    elif parser.getName() == 'save_run':
        try:
            save_run(cmd_pointer, parser)
        except Exception as e:
            return False
        return True
    elif parser.getName() == 'list_runs':
        return list_runs(cmd_pointer, parser)
    elif parser.getName() == 'display_run':
        return display_run(cmd_pointer, parser)
    elif parser.getName() == 'exec_run':
        execute_run(cmd_pointer, parser)

    # File system commands
    elif parser.getName() == 'list_files':
    
        return list_files(cmd_pointer, parser)
    elif parser.getName() == 'import_file':
        return import_file(cmd_pointer, parser)
    elif parser.getName() == 'export_file':
        return export_file(cmd_pointer, parser)
    elif parser.getName() == 'copy_file':
        return copy_file(cmd_pointer, parser)
    elif parser.getName() == 'remove_file':
        return remove_file(cmd_pointer, parser)

    # General commands
    elif parser.getName() == 'welcome':
        # Triggered by `adccl`
        # For testing
        # print(splash(raw=True))
        # print('- - - - - - - - - - - - - -')
        # print(tags_to_markdown(splash(raw=True)))
        # print('- - - - - - - - - - - - - -')
        return output_text(splash(raw=True), nowrap=True)
    elif parser.getName() == 'get_status':
        return return_context(cmd_pointer, parser)
    elif parser.getName() == 'display_history':  # Addressed
        return display_history(cmd_pointer, parser)
    elif parser.getName() == 'display_data':
        return display_data(cmd_pointer, parser)
    elif parser.getName() == 'clear_sessions':
        return clear_other_sessions(cmd_pointer, parser)
    elif parser.getName() == 'edit_config':
        return edit_card(cmd_pointer, parser)

    # Help commands
    elif parser.getName() == 'intro':
        # For testing
        # print(adccl_intro)
        # print('- - - - - - - - - - - - - -')
        # print(tags_to_markdown(adccl_intro))
        # print('- - - - - - - - - - - - - -')
        return output_text(adccl_intro)
    elif parser.getName() == 'docs':
        return docs(cmd_pointer, parser)

    # Show molecules commands
    elif parser.getName() == 'show_molecules':
        return display_mols(cmd_pointer, parser)
    elif parser.getName() == 'show_api_molecules':
        return display_mols(cmd_pointer, parser)

    # Toolkit execution
    elif str(parser.getName()).startswith('toolkit_exec_'):
        try:
            return execute_tookit(cmd_pointer, parser)
        except BaseException as err:
            err = err + '\n' + str(parser.asList())
            return output_error(msg('fail_toolkit_exec_cmd'), cmd_pointer)
    
    return


# Initialises the metadata and workspace directorys for the tool when first run
# if a directory has been deleted it will recreate it
def initialise():
    if not os.path.isdir(_meta_dir):
        os.mkdir(_meta_dir)
    if not os.path.isdir(_meta_dir_toolkits):
        os.mkdir(_meta_dir_toolkits)
    if not os.path.isdir(_meta_workspaces + '/DEFAULT'):
        os.mkdir(_meta_workspaces)
        os.mkdir(_meta_workspaces + '/DEFAULT')
    if not os.path.isfile(_meta_registry):
        initialise_registry()
    if not os.path.isdir(os.path.dirname(_meta_registry_session)):
        os.mkdir(os.path.dirname(_meta_registry_session))
    if not os.path.isfile(_meta_login_registry):
        login_manager.initialise_toolkit_login()


# Intialize the stateful pickle for the user's sessions.


# Open documentation webpage.
def docs(cmd_pointer, parser):
    import webbrowser
    # url = 'data:text/html,<html style="height:100%"><body contenteditable style="height:100%;font-family:sans-serif;font-size:13px;color:dimgray;display:flex;align-items:center;justify-content:center">This is a placeholder for the adccl documentation.</body></html>'
    url = 'https://research.ibm.com/topics/accelerated-discovery'
    webbrowser.open_new(url)
    return output_warning(
        'Our documentation website is yet to be built,\nbut this command will open it in the browser.',
        cmd_pointer,
    )

# adds a registry Toolkit Directory
# in future it will take a tar file and explode it into the directory...
# future work on packaging etc...


def welcome(cmd_pointer, parser):
    """ Display welcome screen """
    return output_text(splash(), nowrap=True)


def return_context(cmd_pointer, parser):
    status = msg('status', cmd_pointer.settings['workspace'], str(cmd_pointer.settings['context']))
    return output_text(status, cmd_pointer, nowrap=True, pad=1)


# List the installed toolkits.
def list_toolkits(cmd_pointer, parser):
    toolkits = []
    table_headers = ('Toolkit', 'Description')

    # Assemble table data.
    for name in cmd_pointer.settings['toolkits']:
        description = load_toolkit_description(cmd_pointer, name)
        toolkits.append(list([name, description]))

    # No toolkits installed yet.
    if len(toolkits) == 0:
        return output_text(msg('no_toolkits_installed'), cmd_pointer, pad=1)

    # Display/return table.
    return output_table(toolkits, cmd_pointer, headers=table_headers, note=msg('all_toolkits_currently_installed'))


# List all available toolkits
def list_all_toolkits(cmd_pointer, parser):
    # This will need to be replaced with a scan of the toolkits directory.
    toolkits = []
    table_headers = ('Toolkit', 'Installed', 'Description')

    # Assemble table data.
    for name in _all_toolkits:
        is_installed = 'Yes' if name in cmd_pointer.settings['toolkits'] else '-'
        description = load_toolkit_description(cmd_pointer, name)
        toolkits.append(list([name, is_installed, description]))

    # Add styling tags
    for i, row in enumerate(toolkits):
        is_installed = row[1] == 'Yes'
        if not is_installed:
            for j, col_text in enumerate(row):
                toolkits[i][j] = f'<soft>{col_text}</soft>'

    # Display/return table.
    return output_table(toolkits, cmd_pointer, headers=table_headers)


# Set the context of the application to and existing toolkit.
# This means the user will only receive access to base commands
# and the toolkit commands of the toolkit currently in context.
def set_context(cmd_pointer, parser):
    if parser['toolkit_name'] is None:
        return return_context(cmd_pointer, parser)
    
    reset = False
    if 'reset' in parser:
        reset=True

    toolkit_name = parser['toolkit_name'].upper()
    toolkit_current = None
    
    if toolkit_name.upper() not in cmd_pointer.settings['toolkits']:
       
        if toolkit_name is None:
            return return_context(cmd_pointer, parser)
        
        # Toolkit doesn't exist.
        return output_error(
            msg('fail_toolkit_not_installed', toolkit_name, split=True),
            cmd_pointer,
            nowrap=True
        )
       
    else:
        
        old_cmd_pointer_context = cmd_pointer.settings['context']
        old_toolkit_current = cmd_pointer.toolkit_current
        load_ok, toolkit_current = load_toolkit(toolkit_name)
        
        if load_ok:
            cmd_pointer.settings['context'] = toolkit_name
            cmd_pointer.toolkit_current = toolkit_current
            refresh_prompt(cmd_pointer.settings)
            write_registry(cmd_pointer.settings, cmd_pointer)
            create_statements(cmd_pointer)
            cmd_pointer.current_help.reset_help()
            #cmd_pointer.current_help.help_current.extend(toolkit_current.methods_help)
            login_success = False
            expiry_datetime = None
            try:
                
                # raise BaseException('Error message') # For testing

                login_success, expiry_datetime = login_manager.load_login_api(cmd_pointer, toolkit_name,reset=reset)

            except BaseException as err:
                # Error logging in.
                output_error(msg('err_login', toolkit_name, err, split=True), cmd_pointer)
                cmd_pointer.settings['context'] = old_cmd_pointer_context
                cmd_pointer.toolkit_current = old_toolkit_current
                return
            
            if not login_success:
                # Failed to log in.
                output_error(msg('err_login', toolkit_name, split=True), cmd_pointer)
                cmd_pointer.settings['context'] = old_cmd_pointer_context
                cmd_pointer.toolkit_current = old_toolkit_current
                return
            
            # Success switching context & loggin in.
            if cmd_pointer.notebook_mode or not cmd_pointer.api_mode:
                return output_success(msg('success_login', toolkit_name, expiry_datetime, split=True), cmd_pointer)
            else:
                return output_text(splash(toolkit_name, cmd_pointer), nowrap=True)

        else:
            # Failed to load the toolkit
            cmd_pointer.settings['context'] = old_cmd_pointer_context
            cmd_pointer.toolkit_current = old_toolkit_current
            return output_error(msg('fail_toolkit_loading', toolkit_name), cmd_pointer)


# Unset the context of the application.
def unset_context(cmd_pointer, parser):
    if cmd_pointer.settings['context'] is None:
        return output_text(msg('no_context_set'), cmd_pointer, pad_btm=1)
    cmd_pointer.settings['context'] = None
    cmd_pointer.toolkit_current = None
    cmd_pointer.current_help.reset_help()
    write_registry(cmd_pointer.settings, cmd_pointer)
    refresh_prompt(cmd_pointer.settings)
    create_statements(cmd_pointer)
    print('')


# Display history of commands.
def display_history(cmd_pointer, parser):
    # readline.write_history_file(cmd_pointer.histfile)  # @Phil, I put this back but it was commented out

    history = []

    # Fetch last 30 commands from history.
    i = 0
    hist_len = readline.get_current_history_length()
    index_col_width = 3
    min_gap = 2
    reached_bottom = False

    while i < 30:
        # Make sure all items are aligned.
        line_nr = hist_len - i
        line_nr_length = len(str(line_nr))
        if line_nr_length + min_gap > index_col_width:
            index_col_width = line_nr_length + min_gap
        gap = (index_col_width - line_nr_length) * ' '

        # Fetch history item.
        try:
            entry = readline.get_history_item(line_nr)
            i = i + 1
            if entry:
                entry = entry.replace('\n', '')
                entry_str = f'<cmd>{entry}</cmd>'
            else:
                entry_str = '<soft>Workspace created</soft>'
            history.append((line_nr, entry_str))
            
            
            # Reached bottom.
            if not entry:
                reached_bottom = True
                break
        except BaseException as err:
            output_error(msg('err_fetch_history', err, split=True), cmd_pointer)
            i = 31

    # history = list(map(lambda t: [f'<soft>{t[0]}</soft>{t[1]}{t[2]}'], history)) # Using table - trash

    # Add ellipsis if history is longer than 30 items.
    if not reached_bottom:
        history.append((line_nr - 1,  '<soft>...</soft>'))
    history.reverse()

    # Display/return table.
    
    return output_table(history, cmd_pointer, headers=['',  'Command History'])

    # Trash
    # if cmd_pointer.notebook_mode == False:
    #     history = list(map(lambda t: f'<soft>{t[0]}</soft>{t[1]}{t[2]}', history))
    #     history = '<h1>Command History</h1>\n' + '\n'.join(history)
    #     return output_text(history, cmd_pointer, pad=1, nowrap=True)

    # else:
    #     import pandas
    #     from helpers.style_parser import strip_tags
    #     pandas.set_option('display.max_colwidth', None)
    #     history = list(map(lambda x: [x[0], x[1], strip_tags(x[2])], history))
    #     return pandas.DataFrame(history, columns=['', '', 'Command History'], index=None)


# Display a csv file in a table.
def display_data(cmd_pointer, parser):
    import pandas
    workspace_path = cmd_pointer.workspace_path(cmd_pointer.settings['workspace'].upper()) + '/'
    file_path = parser['file_path']
    filename = file_path.split('/')[-1]

    # Allow for no extension.
    if len(filename.split('.')) == 1:
        filename = filename + '.csv'
        file_path = file_path + '.csv'

    try:
        if filename.split('.')[-1].lower() == 'csv':
            # From csv file.
            try:
                df = pandas.read_csv(workspace_path + file_path)
                return output_table(df)
            except FileNotFoundError:
                return output_error(msg('fail_file_doesnt_exist', file_path), cmd_pointer)
            except BaseException as err:
                return output_error(msg('err_load_csv', err, split=True), cmd_pointer)
        else:
            return output_error(msg('invalid_file_format', 'csv', split=True), cmd_pointer)

    except BaseException as err:
        output_error(msg('err_unknown', err, split=True), cmd_pointer)


def edit_card(cmd_pointer, parser):
    if cmd_pointer.notebook_mode == True:
        print('Json Editor only available from command line')
        return True
    import  ad4e_opentoolkit.edit_json.edit_json as edit_json
    file_to_edit = cmd_pointer.workspace_path(
        cmd_pointer.settings['workspace'].upper()) + '/' + parser.as_dict()['json_file']
    if not os.path.isfile(file_to_edit):
        print('file does not exist')
        return False
    edit_json(file_to_edit)
    return True
