#!/usr/local/opt/python@3.9/bin/python3.9
# Copyright 2022 IBM, Inc. or its affiliates. All Rights Reserved.

import os
import logging
import sys
import readline
import ad4e_opentoolkit.app.login_manager
import string
import uuid
from cmd import Cmd

# Main
from ad4e_opentoolkit.app.main_lib import lang_parse, initialise, set_context
from ad4e_opentoolkit.toolkit.toolkit_main import load_toolkit

# Core
# from core.grammar import *
import ad4e_opentoolkit.core.help as adccl_help
from ad4e_opentoolkit.core.grammar import grammar_help, statements, statements_def, create_statements,output_train_statements
from ad4e_opentoolkit.core.lang_sessions_and_registry import write_registry, load_registry, delete_session_registry
from ad4e_opentoolkit.core.lang_workspaces import set_workspace
import ad4e_opentoolkit.app.login_manager as login_manager
# Helpers
from ad4e_opentoolkit.helpers.general import singular, confirm_prompt
from ad4e_opentoolkit.helpers.output import msg, output_text, output_error, output_warning
from ad4e_opentoolkit.helpers.general import refresh_prompt
from ad4e_opentoolkit.helpers.splash import splash

# Globals
from ad4e_opentoolkit.app.global_var_lib import _repo_dir as _repo_dir
from ad4e_opentoolkit.app.global_var_lib import _meta_dir as _meta_dir
from ad4e_opentoolkit.app.global_var_lib import _meta_workspaces as _meta_workspaces
from ad4e_opentoolkit.app.global_var_lib import _meta_login_registry as _meta_login_registry
from ad4e_opentoolkit.app.global_var_lib import _all_toolkits as _all_toolkits
from ad4e_opentoolkit.app.global_var_lib import  _meta_dir_toolkits as _meta_dir_toolkits



sys.ps1 = '\x01\033[31m\x02>>> \x01\033[0m\x02'


# Used for for converting lists to strings.
def convert(lst):
    return str(lst).translate('[],\'')


# this is the command class/object that is the center of the command line DSL Shell environment
# it holds the Parsed grammar and current state of a users engagement in the utility
class run_cmd(Cmd):
    space = " "
    # __all__                     =   ["Cmd"]
    # PROMPT                      =   '(Cmd) '
    IDENTCHARS              = string.ascii_letters + string.digits + '_'
    intro                   = '/'  # This is defined in cmdloop() below.
    home_dir                =_meta_dir
    repo_dir                =_repo_dir
    current_statements      = statements
    current_statement_defs  = statements_def
    toolkit_dir             =_meta_dir_toolkits
    complete_index          = None
    complete_orig_line      = None
    
    settings                = None
    original_settings       = None
    session_id              = '_session_' + str(uuid.uuid4()).replace("-", "")
    toolkit_current         = None
    prompt                  = None
    histfile                = os.path.expanduser(_meta_dir + '/.cmd_history')
    histfile_size           = 1000
    current_help            = adccl_help.adccl_help()
    current_help.help_orig  = grammar_help.copy()
    current_help.reset_help()
    notebook_mode           = False
    api_mode                = False
    login_settings          = None
    api_variables           = {}
    llm_handle              = None
    refresh_vector          = False
    refresh_train           = False
    llm_service             = 'OPENAI'
    llm_model               = 'gpt-3.5-turbo'
    llm_models              =   {'OPENAI':'gpt-3.5-turbo','WATSONX':'mosaicml/mpt-7b'}

    def workspace_path(self, workspace: str):
        try:
            x = os.path.expanduser(self.settings['paths'][workspace.upper()] + '/' + workspace.upper())

            return x
        except BaseException:
            return os.path.expanduser(_meta_workspaces + '/' + workspace.upper())

    def set_workspace_path(self, workspace: str, path: str):
        self.settings['paths'][workspace.upper()] = str

    # Initialises the Class for Run command
    def __init__(self, completekey='Tab',notebook=False,api=False):
        self.notebook_mode=notebook
        self.api_mode=api
        super().__init__()
        if sys.platform == 'darwin':
            if 'libedit' in readline.__doc__:
                readline.parse_and_bind("bind ^I rl_complete")
            else:
                readline.parse_and_bind("tab: complete")
        readline.set_completer(self.complete)
        
        self.settings = load_registry(self, orig_reg=True)
        self.original_settings = load_registry(self, orig_reg=True)
        write_registry(self.settings, self)
        
        self.prompt = refresh_prompt(self.settings)
        
        if self.settings['context'] in self.settings['toolkits']:
            ok, toolkit_current = load_toolkit(self.settings['context'])
            if ok:
                self.toolkit_current = toolkit_current
                create_statements(self)
        
        self.login_settings = login_manager.load_login_registry()
        
        if self.settings['workspace'] is not None:
            self.histfile = os.path.expanduser(self.workspace_path(self.settings['workspace'].upper()) + '/.cmd_history')
        
        if self.settings['context'] is not None:
            login_manager.load_login_api(self, self.settings['context'])
        
        try:
          if self.settings['env_vars']['refresh_help_ai']==True:
            self.refresh_vector=True
            self.refresh_train=True
        except:
            pass
        try:
            self.llm_service= self.settings['env_vars']['llm_service']
            self.llm_model=self.llm_models[self.llm_service]
                  
        except Exception as e:
            #print(e)
            #print("failed to load service llm")
            pass
        
        output_train_statements(self)



    def do_help(self, inp):
        """
        Displaying help.

        The different entry points:
            ? list                   --> The questionmark is interpreted and stripped by the language parser
            list ?                   --> This goes via the run_cmd.default() function
            python3 main.py '? list' --> This goes via __main__
            %adccl ? list            --> Notebook and API requests go via api_remote()

        """

        # `??` --> Advanced help.
        starts_with =False
        if inp.strip() == '?':
            return output_text(adccl_help.advanced_help(), self, pad=1)

        # Strip question marks at the beginning and end.
        if len(inp.strip()) > 0 and inp.split()[0] == '?':
            # Note: Usually the language parser will strip the question mark,
            # but this is still needed you run ```python3 main.py '? list'```.
            inp = inp.lstrip('?')
        elif len(inp.strip()) > 0 and inp.split()[-1] == '?':
            inp = inp.rstrip('?')
            starts_with=True
        
        inp = inp.lower().strip()
        one_word_cmd = len(inp.split()) == 1
        all_commands = self.current_help.help_current
        
        matching_commands = []

        # `?` --> Display all commands.
        if len(inp.split()) == 0:
            return output_text(
                adccl_help.all_commands(all_commands, toolkit_current=self.toolkit_current, cmd_pointer=self),
                self,
                pad=2,
                tabs=1,
                nowrap=True
            )

        # `<toolkit_name> ?` --> Display all toolkkit commands.
        if inp.upper() in _all_toolkits:
            toolkit_name = inp.upper()
            ok, toolkit = load_toolkit(toolkit_name)
            return output_text(
                adccl_help.all_commands(toolkit.methods_help, toolkit_name, cmd_pointer=self),
                self,
                pad=2,
                tabs=1,
                nowrap=True
            )
        try:
            for i in self.toolkit_current.methods_help:
                if i not in all_commands:
                    all_commands.append(i)
        except:
            pass
      
        # Look for commands that include this exact word singular/plural.
        if one_word_cmd and starts_with==False:
            for command in all_commands:
                command_str = str(command['command']).strip().lower()
                words = command_str.split()
                inp_singular = singular(inp)
                inp_plural = inp_singular + 's'
                if inp_singular in words or inp_plural in words:
                    query_type = 'word_match'
                    matching_commands.append(command)

        # Look for commands that start with this string.
        if not len(matching_commands):
            
            for command in all_commands:
                command_str = str(command['command']).strip().lower()
                if command_str.upper().startswith(inp.upper()):
                    query_type = 'starts_with'
                    matching_commands.append(command)

        if len(matching_commands) > 1:
            return output_text(adccl_help.queried_commands(matching_commands, inp=inp, query_type=query_type), self, pad=1, nowrap=True)
        elif len(matching_commands) == 1:
            return output_text(adccl_help.command_details(matching_commands[0]), self, pad=2, tabs=1)
        else:
            return output_error(msg('err_invalid_cmd', None, split=True), self)

    # Preloop is called by cmd to get an update the history file
    # Each History File
    def preloop(self):
        if readline and os.path.exists(self.histfile):
            try:
                readline.read_history_file(self.histfile)
            except BaseException:
                # Create history file in case it doesn't exist yet.
                # - - -
                # To trigger:
                # >> create new workspace foobar
                # >> ctrl+c
                # (Reboot)
                readline.write_history_file(self.histfile)

    # Post loop is called by cmd to get an update the history file
    def postloop(self):
        if readline:
            readline.set_history_length(self.histfile_size)
            readline.write_history_file(self.histfile)

    def add_history(self, inp):
        #readline.add_history('\001'+inp+'\002')
        readline.add_history(inp)

    #####################################################################################################################
    # This is the Auto Complete Method that gets called on Forward Tab.
    # Currently parsing of the pyparsing statements and finding the statement that it fails against at a character furtheresrt long the statement string is the method used.
    # this is an area for improvement further along the line.
    def match_display_hook(self, substitution, matches, longest_match_length):
        print('\n----------------------------------------------\n')
        for match in matches:
            print(match)
        print(self.prompt.rstrip() + readline.get_line_buffer())
        readline.redisplay()

    def complete(self, text, state):
        if state == 0:
            orig_line = readline.get_line_buffer()

        started_command = None
        i_s = 0
        yy = []

        if len(orig_line.split()) > 1:
            orig_word = orig_line.split()[len(orig_line.split()) - 1]
        else:
            orig_word = orig_line
        # print(orig_line)
        # print(self.current_statements.CloseMatch(orig_line))
        matches = []
        # set_completion_display_matches_hook(self.match_display_hook)  #Look at later
        test_list = []

        while len(yy) == 0 and i_s < len(self.current_statements):

            a, b = self.current_statements[i_s].run_tests(orig_line, printResults=False, fullDump=False)
            test_list.append(b[0])
            i_s = i_s + 1
        best_fit = 0

        for x in test_list:
            if error_col_grabber(str(x)) > best_fit:
                best_fit = error_col_grabber(str(x))

        for i in test_list:
            if error_col_grabber(str(i)) < best_fit:
                continue

            if len(i) > 1:
                c = i[1]
                x = c.explain()

                x = x.replace(orig_line, '')
                if (x.split(',')[0].find("Expected CaselessKeyword") > -1 or x.split(',')
                        [0].find("Expected Keyword") > -1) and x.split(',')[0].find("'" + orig_word.lower()) > -1:

                    yy = x.split(",")[0].split("'")[1]
                    readline.insert_text(yy[len(orig_word):])
                    #readline.redisplay()

                    readline.insert_text(" ")
                    readline.redisplay()

                    return '' 
        for i in test_list:
            if error_col_grabber(str(i)) < best_fit:
                continue

            if len(i) > 1:
                c = i[1]
                x = c.explain()
                x = x.replace(orig_line, '')
                if (x.split(',')[0].find("Expected CaselessKeyword") > -1 or x.split(',')
                        [0].find("Expected Keyword") > -1) and x.split(',')[1].find("at char 0") > -1:

                    if str(str(i[1]).split(",")[0].split("Keyword")[1].split("'")[1]).strip().upper() == str(
                            i[0] + x.split(',')[0].split("Keyword")[1].split("'")[1]).strip().upper():
                        
                        readline.insert_text(x.split(',')[0].split("Keyword")[1].split("'")[1].strip())
                        #readline.redisplay()

                        readline.insert_text(" ")
                        readline.redisplay()

                        return ''
                    continue

                # print(x)

        for i in test_list:
            if error_col_grabber(str(i)) < best_fit:
                continue

            if len(i) > 1:
                c = i[1]
                x = c.explain()
                x = x.replace(orig_line, '')
                if (x.split(',')[0].find("Expected CaselessKeyword") > -1 or x.split(',')
                        [0].find("Expected Keyword") > -1) and x.split(',')[0].find("'" + orig_word.lower()) == -1:
                    spacing = ""
                    if len(orig_line) == len(i[0]):
                        spacing = " "
                   
                    if error_col_grabber(x)-1 < len(orig_line):
                        if  len(orig_line[error_col_grabber(x)-1:len(orig_line)].strip())> 0:
                            return []
                    readline.insert_text(spacing + x.split(',')[0].split("Keyword")[1].split("'")[1].strip())
                    #readline.redisplay()

                    readline.insert_text(" ")
                    readline.redisplay()

                    return ''

        for i in test_list:
            if error_col_grabber(str(i)) < best_fit:
                continue

            if len(i) > 1:
                c = i[1]
                x = c.explain()
                x = x.replace(orig_line, '')

                if (x.split(',')[0].find("Expected '('") > -1 or x.split(',')[0].find("Expected ')'") > -1):
                    # print('here')
                    # print(x)
                    if x.split(',')[0].find("Expected '('") > -1:
                        readline.insert_text('(')
                    else:
                        readline.insert_text(')')

                    readline.redisplay()

                    return ''
                if x.split(',')[0].find("Expected string enclosed in '\"'"):
                    # print('here')
                    # print(x)

                    readline.insert_text("'")
                    readline.redisplay()

                    return ''
                else:
                    pass
                    # xx=x[x.find("{")+1:x.find('}')].split('|')

        return []

    # Catches the exit command

    def do_exit(self,dummy_inp_do_not_remove):
        write_registry(self.settings, self, True)
        delete_session_registry(self.session_id)
        # readline.remove_history_item(readline.get_current_history_length()-1)
        '''exiting the application. Shorthand: x q.'''
        return True

    # prevents on return on a blank line default receiving the previous input to run

    def emptyline(self):
        pass

    # Default method call on hitting of the return Key, it tries to parse and execute the statements.
    def default(self, inp):
        x = None
       
        if convert(inp).split()[-1] == '?' and not convert(inp).upper().startswith('TELL ME'):
            
            return self.do_help(inp)
        
        try:
            try:
                self.settings = load_registry(self)
            except:
                # Brutal situation where someone hit clear sessions in another session , shut down abruptly so as not to kill registry file
                print('fatal error session registry not avaiable, performing emergency shutdown !')
                self.do_exit('exit emergency')
            
            y = self.current_statement_defs.parseString(convert(inp), parseAll=True)
            
            x = lang_parse(self, y)
        
            self.prompt = refresh_prompt(self.settings)
            logging.info('Ran: ' + inp)
        except BaseException as err1:
            #print(err1)
            # Removing due to usability being able to recall item and correct:
            # try:
            #    readline.remove_history_item(readline.get_current_history_length()-1) # Does not save an incorrect instruction
            # except:
            #    Error_descriptor=None
            error_descriptor = None
            error_col = -1
            invalid_command = False
            i_s = 0
            while i_s < len(self.current_statements):
                a, b = self.current_statements[i_s].runTests(convert(inp), printResults=False, fullDump=False, parseAll=True)

                for i in b:
                    if len(i) > 1:
                        invalid_command = True

                        # Fetch error description
                        # @Phil, these are not very helpful, is the idea here to make suggestions for longer almost-correct commands?
                        c = i[1]
                        try:
                            x = c.explain()
                        except:
                            print("unknown Error: ")
                            print(err1)
                        

                        if x.find("Expected CaselessKeyword") > -1 and x.find('at char 0') == -1:
                            if error_col < error_col_grabber(x):
                                error_descriptor = x.replace("CaselessKeyword", "keyword").replace("ParseException:", "Syntax Error:: ")
                                error_col = error_col_grabber(x)
                                
                        elif x.find("found end of text") > -1 and x.find('at char 0') == -1:
                            if error_col < error_col_grabber(x):
                                error_descriptor = x.replace("ParseException:", "Syntax Error:: ")
                                error_col = error_col_grabber(x)
                        # @Phil these general errors tend to be unhelpful.
                        #@moenen, on back log to improve but for a lot of users familiar from database command line errors
                        else:
                            if error_col < error_col_grabber(x):
                                 error_descriptor = x
                                 error_col = error_col_grabber(x)
                i_s = i_s + 1

            # Print error
            if invalid_command:
                # Note: error_descriptor is optional.
                if error_col_grabber(error_descriptor) == 0:
                    if self.notebook_mode is True:
                        return output_error(msg('err_invalid_cmd', 'Not a Valid Command, try "?" to list valid commands', split=True), self)
                    else:
                        output_error(msg('err_invalid_cmd', 'Not a Valid Command, try "?" to list valid commands', split=True), self)
                else:
                    error_msg=error_descriptor.split("Syntax")[0]
                    if self.notebook_mode is True:
                        from IPython.display import display
                        if error_col_grabber(error_descriptor) == 1:
                            display(output_error(msg('err_invalid_cmd', error_msg.split('Parse')[0], split=True),return_val=True, cmd_pointer=self))
                            display(output_text("Perhaps you could try one of the following:",return_val=True,cmd_pointer=self))
                            display(self.do_help(error_first_word_grabber(error_descriptor)+' ?'))
                        else:
                            display(output_error(msg('err_invalid_cmd', error_msg, split=True),return_val=True, cmd_pointer=self))
                            display(output_text("Perhaps you could try one of the following:",return_val=True,cmd_pointer=self))
                            display(self.do_help(inp[0:error_col_grabber(error_descriptor)-1]+' ?'))
                        return output_text("If there is not an option that meets your requirement type '?' to list all command options",return_val=True,cmd_pointer=self)
                        
                    else:
                        
                        
                        if error_col_grabber(error_descriptor) == 1:
                            output_error(msg('err_invalid_cmd', error_msg.split('Parse')[0], split=True), self)
                            output_text("Perhaps you could try one of the following:",self)
                            self.do_help(error_first_word_grabber(error_descriptor)+' ?')
                        else:
                            output_error(msg('err_invalid_cmd', error_msg, split=True), self)
                            output_text("Perhaps you could try one of the following:",self)                     
                            self.do_help(inp[0:error_col_grabber(error_descriptor)-1]+' ?')
                        
                        output_text("If there is not an option that meets your requirement type '?' to list all command options",return_val=True,cmd_pointer=self)
                return False


            else:
                output_error(msg('err_invalid_cmd', x, split=True), self)
                return False
                #return output_error(msg('err_unknown', err1, split=True), self) # @moenen this was not catching the error returned by the function and re-issuing splash screen
        if self.refresh_train ==True:
            output_train_statements(self)
            self.refresh_train=False
        if self.notebook_mode is True:
            return x
        elif self.api_mode==False:
            if x not in (True,False,None):    
                
                print(x)
            else:
                return


# this function retuns the errror positioning in the statement tht has been parsed.
def error_col_grabber(error):
    e = error.split('col:')[1]
    e1 = e.replace(')', '')
    return int(e1)

def error_first_word_grabber(error):
    word = error.split('found ')[1].split("'")[1]
   
    return str(word)


# Main execution application
# if the application is called with parameters it executes as Parameters,
# if called without parameters the command line enters the shell environment
# History is only kept for commands executed once in the shell
def api_remote(inp: str, connection_cache: dict = _meta_login_registry, api_context: dict = {'workspace': None, 'toolkit': None}, api_var_list={}):
    initialise()
    arguments = inp.split()
    inp = ''
    a_space = ''
    
    magic_prompt = run_cmd(notebook=True)
    
    connection_cache = magic_prompt.login_settings
    magic_prompt.notebook_mode = True
    if api_context['workspace'] is None:
        api_context['workspace'] = magic_prompt.settings['workspace']
    else:
        x = {'Workspace_Name': api_context['workspace']}
        set_workspace(magic_prompt, x)

    if api_context['toolkit'] is None:
        api_context['toolkit'] = magic_prompt.settings['context']
    else:
        x = {'toolkit_name': api_context['toolkit']}
        set_context(magic_prompt, x)
        
    magic_prompt.api_variables = api_var_list
    try:
        readline.read_history_file(magic_prompt.histfile)
    except BaseException:
        readline.add_history('')
        readline.write_history_file(magic_prompt.histfile)
        readline.read_history_file(magic_prompt.histfile)
    for i in arguments:
        inp = inp + a_space + i
        a_space = ' '

    # Check for Help from a command line request
    if len(inp.strip()) > 0:
        if inp.split()[0] == '?' or (inp.split()[-1] == '?' and not convert(inp).upper().startswith('TELL ME')) or inp.strip() == '??':
            if inp.strip() == '?':
                inp = ''
            elif inp.strip() == '??':
                inp = '?'

            return magic_prompt.do_help(inp)
        else:
            # if there is a argument and it is not help attemt to run the command
            # Note, may be possible add code completion here #revisit

            magic_prompt.preloop()
            magic_prompt.add_history(inp)
            magic_prompt.postloop()
            readline.write_history_file(magic_prompt.histfile)

            result = magic_prompt.default(inp)
            api_context['workspace'] = magic_prompt.settings['workspace']
            api_context['toolkit'] = magic_prompt.settings['context']

            magic_prompt.do_exit('dummy do not remove')

            return result

def cmd_line():
    initialise()
    inp = ''
    a_space = ''

    for i in sys.argv[1:]:
        inp = inp + a_space + i
        a_space = ' '

    command_line = run_cmd()
    # Check for help from a command line request.
    if len(inp.strip()) > 0:
        words = inp.split()
        if inp.split()[0] == '?' or inp.split()[-1] == '?' or inp.strip() == '??':
            if inp.strip() == '?':
                inp = ''
            elif inp.strip() == '??':
                inp = '?'
        # if words[0] == '?' or words[-1] == '?':
            command_line.do_help(inp.strip())
        elif words[0].lower() == '-s':
            for i in words:
                
                print(i)
            set_workspace(command_line, {'Workspace_Name': words[1].upper()})
            set_context(command_line, {'toolkit_name': words[2].upper()})
            command_line.preloop()
            command_line.add_history(str(' '.join(words[3:])).strip())
            command_line.postloop()
            command_line.default(str(' '.join(words[3:])).strip())
        else:
            # if there is a argument and it is not help attemt to run the command
            # Note, may be possible add code completion here #revisit
            command_line.preloop()
            command_line.add_history(inp.strip())
            command_line.postloop()
            command_line.default(inp.strip())
        command_line.do_exit('dummy do not remove')
    else:
        # If no argument passed then enter
        exit = False
        while exit == False:
            try:
                # The cmdloop parameter controls the startup screen, it overrides self.intro.
                command_line.cmdloop(splash(command_line.settings['context'], command_line, startup=True))
                exit = True
            except KeyboardInterrupt:
                command_line.postloop()
                if confirm_prompt('Are you sure you wish to exit?'):
                    exit = True
                    command_line.do_exit('dummy do not remove')
            except BaseException as err:
                output_error(msg('err_invalid_cmd', err, split=True), command_line)

    
if __name__ == "__main__":
    cmd_line()
    