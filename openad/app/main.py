"""Main application contain runtime class RUNCMD()"""
#!/usr/local/opt/python@3.9/bin/python3.9
# Copyright 2022 IBM, Inc. or its affiliates. All Rights Reserved.

import os
import sys
import readline
import re
import string
import uuid
from cmd import Cmd

# Main
from openad.app.main_lib import lang_parse, initialise, set_context, unset_context
from openad.toolkit.toolkit_main import load_toolkit
from openad.app import login_manager

# Core
import openad.core.help as openad_help
from openad.core.grammar import (
    grammar_help,
    statements,
    statements_def,
    create_statements,
    output_train_statements,
)
from openad.core.lang_sessions_and_registry import write_registry, load_registry, delete_session_registry
from openad.core.lang_workspaces import set_workspace
from openad.app.memory import Memory


from openad.llm_assist.model_reference import SUPPORTED_TELL_ME_MODELS, SUPPORTED_TELL_ME_MODELS_SETTINGS

# Helpers
from openad.helpers.general import singular, confirm_prompt
from openad.helpers.output import msg, output_text, output_error, output_warning, strip_tags
from openad.helpers.general import refresh_prompt
from openad.helpers.splash import splash
from openad.helpers.output_content import info_workspaces, info_toolkits, info_runs, info_context

# Globals
from openad.app.global_var_lib import _repo_dir
from openad.app.global_var_lib import _meta_dir
from openad.app.global_var_lib import _meta_workspaces
from openad.app.global_var_lib import _all_toolkits
from openad.app.global_var_lib import _meta_dir_toolkits


sys.ps1 = "\x01\033[31m\x02>>> \x01\033[0m\x02"


# Used for for converting lists to strings.
def convert(lst):
    """Used for for converting lists to strings."""
    return str(lst).translate("[],'")


class RUNCMD(Cmd):
    """
    The center of the command line DSL Shell environment, this holds the parsed
    grammar and current state of a user's engagement.
    """

    space = " "
    IDENTCHARS = string.ascii_letters + string.digits + "_"
    intro = "/"  # This is defined in cmdloop() below.
    home_dir = _meta_dir
    repo_dir = _repo_dir
    current_statements = statements
    current_statement_defs = statements_def
    toolkit_dir = _meta_dir_toolkits
    complete_index = None
    complete_orig_line = None
    settings = None
    original_settings = None
    session_id = "_session_" + str(uuid.uuid4()).replace("-", "")
    toolkit_current = None
    prompt = None
    histfile = os.path.expanduser(_meta_dir + "/.cmd_history")
    histfile_size = 1000  # prompt history file per workspace limit
    current_help = openad_help.OpenadHelp()  # handle to the current help object
    current_help.help_orig = grammar_help.copy()  # copy of the base line command help functions (excludes Toolkits)
    current_help.reset_help()  # initialises help
    notebook_mode = False  # set to denote if the calls are coming from a jupyter notebook
    api_mode = False  # set to denote app is called from an external API
    login_settings = None  # where the login settings get intitialised to
    api_variables = {}  # variables passed to from external applications
    # Servicing the LLM related Function States
    llm_handle = None  # connection handle for LLM for Tell Me and other functions
    refresh_vector = (
        False  # Signals the Refresh of the vector DB should be done due to changes in Workspace or Toolkits
    )
    refresh_train = False  # Signals Refreshing of the training repository for help should be done
    llm_service = "OPENAI"  # set with OPENAI as default type until WatsonX or alternative available
    llm_model = "gpt-3.5-turbo"
    llm_models = SUPPORTED_TELL_ME_MODELS_SETTINGS
    # Instantiate memory class
    memory = Memory()

    # Instantiate list of Molecuels for reference
    molecule_list = []
    last_external_molecule = None

    def workspace_path(self, workspace: str):
        """Returns the default workspace directory path"""
        try:
            x = os.path.expanduser(self.settings["paths"][workspace.upper()] + "/" + workspace.upper())
            return x
        except Exception:  # pylint: disable=broad-exception-caught
            # various exceptions can cause this... Any error results in same outcome
            return os.path.expanduser(_meta_workspaces + "/" + workspace.upper())

    def set_workspace_path(self, workspace: str, path: str):
        """Sets the current workspace path in the settings dictionary"""
        self.settings["paths"][workspace.upper()] = os.path.expanduser(path)

    # Initialises the class for Run command.
    def __init__(self, completekey="Tab", notebook=False, api=False):
        self.notebook_mode = notebook
        self.api_mode = api
        super().__init__()

        # This is necessary to ensure readline works predicably and compatibly across MacOS and Linux
        if sys.platform == "darwin":
            if "libedit" in readline.__doc__:
                readline.parse_and_bind("bind ^I rl_complete")
            else:
                readline.parse_and_bind("tab: complete")
        readline.set_completer(self.complete)

        self.settings = load_registry(self, orig_reg=True)  # loads up the sessions settings
        self.original_settings = load_registry(
            self, orig_reg=True
        )  # for reference keeps a copy of original settings on startup

        write_registry(self.settings, self)  # writes the session registry settings

        self.prompt = refresh_prompt(self.settings)  # sets the command prompt

        # load the toolkit in current context
        if self.settings["context"] in self.settings["toolkits"]:
            ok, toolkit_current = load_toolkit(self.settings["context"])
            if ok:
                self.toolkit_current = toolkit_current
                create_statements(self)
        # Initialise current toolkit registry
        self.login_settings = login_manager.load_login_registry()

        # Check for reset in workspace and if so set to default
        if self.settings["workspace"] is not None:
            self.histfile = os.path.expanduser(
                self.workspace_path(self.settings["workspace"].upper()) + "/.cmd_history"
            )
        # set current context login
        if self.settings["context"] is not None:
            success, expiry = login_manager.load_login_api(self, self.settings["context"])
            if success is False:
                self.settings["context"] = None
                self.toolkit_current = None
                unset_context(self, None)
                self.prompt = refresh_prompt(self.settings)
                output_text("Unable to set context on Login, defaulting to no context set.", self, return_val=False)
        try:
            if self.settings["env_vars"]["refresh_help_ai"] is True:
                self.refresh_vector = True
                self.refresh_train = True
        except Exception:  # pylint: disable=broad-exception-caught  # if LLM not initiated move on
            pass
        # Try to load variables for llm. If missing, just pass and move on.
        try:
            self.llm_service = self.settings["env_vars"]["llm_service"]
            self.llm_model = self.llm_models[self.llm_service]
        except Exception:  # pylint: disable=broad-exception-caught  # if LLM not initiated move on
            pass

        output_train_statements(self)

    def do_help(self, inp, display_info=True, starts_with_only=False, **kwargs):
        """CMD class called function:
        Display help about a command, for example 'list'.

        Parameters:
            inp: The input string.
            display_info: If True, display info text when relevant.
                Eg. `workspaces ?`, `toolkits ?`, `runs ?`, etc.

        The different entry points:
            ? list                   --> The questionmark is interpreted and stripped by the language parser
            list ?                   --> This goes via the RUNCMD.default() function
            python3 main.py '? list' --> This goes via __main__
            %openad ? list           --> Notebook and API requests go via api_remote()
        """

        # `??` --> Advanced help (to be implemented)
        if inp.strip() == "?":
            return output_text(openad_help.advanced_help(), self, pad=1)

        # Strip question marks at the beginning and end of input.
        if len(inp.strip()) > 0 and inp.split()[0] == "?":
            # Beginning
            # - - -
            # Note: Usually the language parser will strip the question mark,
            # but this is still needed you run ```python3 main.py '? list'```.
            inp = inp.lstrip("?")
        elif len(inp.strip()) > 0 and inp.split()[-1] == "?":
            # End
            starts_with_only = True
            inp = inp.rstrip("?")

        inp = inp.lower().strip()
        # [:] is to make a copy of the list, so we don't modify the original.
        all_commands = self.current_help.help_current[:]
        matching_commands = {
            "match_word": [],
            "match_start": [],
            "match_anywhere": [],
        }

        # When returning result for a help query, we have to remove
        # the `... ?` and `? ...` commands because their documentation
        # is included in the command string and we don't want them to
        # show up when you query eg. `all ?`. When you just run `?`
        # the inp will be "" and we don't want to remove anything.
        # See grammar.py -> "command help 1" and "command help 2".
        if len(inp):
            cmds_to_ignore = [cmd for cmd in all_commands if "-->" in cmd["command"]]
            for cmd in cmds_to_ignore:
                all_commands.remove(cmd)

        # `?` --> Display all commands.
        if len(inp.split()) == 0:
            return output_text(
                openad_help.all_commands(all_commands, toolkit_current=self.toolkit_current, cmd_pointer=self),
                self,
                pad=2,
                tabs=1,
            )

        # Display info text about important key concepts.
        if display_info and ("return_val" not in kwargs or not kwargs["return_val"]):
            if inp.lower() == "workspace" or inp.lower() == "workspaces":
                output_text("<h1>About Workspaces</h1>\n" + info_workspaces, self, edge=True, pad=1, return_val=False)
            elif inp.lower() == "toolkit" or inp.lower() == "toolkits":
                output_text("<h1>About Toolkits</h1>\n" + info_toolkits, self, edge=True, pad=1, return_val=False)
            elif inp.lower() == "run" or inp.lower() == "runs":
                output_text("<h1>About Runs</h1>\n" + info_runs, self, edge=True, pad=1, return_val=False)
            elif inp.lower() == "context" or inp.lower() == "contexts":
                output_text("<h1>About Context</h1>\n" + info_context, self, edge=True, pad=1, return_val=False)

        # `<toolkit_name> ?` --> Display all toolkkit commands.
        if inp.upper() in _all_toolkits + ["DEMO"]:  # DEMO is omitted from _all_toolkits
            toolkit_name = inp.upper()
            ok, toolkit = load_toolkit(toolkit_name)
            return output_text(
                openad_help.all_commands(toolkit.methods_help, toolkit_name, cmd_pointer=self), self, pad=2, tabs=1
            )

        # Add the current toolkit's commands to the list of all commands.
        try:
            for i in self.toolkit_current.methods_help:
                if i not in all_commands:
                    all_commands.append(i)
        except Exception:  # pylint: disable=broad-exception-caught # do not need to know exception
            pass

        # First list commands with full word matches.
        if not starts_with_only:
            for command in all_commands:
                words = command["command"].split()
                inp_singular = singular(inp)
                inp_plural = inp_singular + "s"
                if inp_singular in words or inp_plural in words:
                    matching_commands["match_word"].append(command)

        # Then list commands starting with the input string.
        for command in all_commands:
            if re.match(re.escape(inp), command["command"]) and command not in matching_commands["match_word"]:
                matching_commands["match_start"].append(command)

        # Then list commands containing the input string.
        if not starts_with_only:
            for command in all_commands:
                if (
                    re.search(re.escape(inp), command["command"])
                    and command not in matching_commands["match_word"]
                    and command not in matching_commands["match_start"]
                ):
                    matching_commands["match_anywhere"].append(command)

        # else:
        all_matching_commands = (
            matching_commands["match_word"] + matching_commands["match_start"] + matching_commands["match_anywhere"]
        )

        result_count = len(all_matching_commands)

        # Check if there is an exact match.
        # This is for case like `run <run_name> ?` which would otherwise
        # display `run <run_name>` as well as `display run <run_name>`
        all_matching_commands_str = [item["command"] for item in all_matching_commands]
        exact_match = inp in all_matching_commands_str

        # This lets us pass a custom padding value.
        # This is used when suggesting commands after your input was not recognized.
        # After "You may want to try:" we don't want a linebreak
        if "pad_top" in kwargs:
            pad = None
            pad_top = kwargs["pad_top"]
            del kwargs["pad_top"]
        else:
            pad = 1
            pad_top = None

        # DISPLAY:
        # No matching commands -> error.
        if result_count == 0 and not exact_match:
            if starts_with_only:
                return output_error(msg("err_no_cmds_starting", inp), self, **kwargs)
            else:
                return output_error(msg("err_no_cmds_matching", inp), self, **kwargs)

        # Single command -> show details.
        elif result_count == 1 or exact_match:
            return output_text(
                openad_help.command_details(all_matching_commands[0], self),
                self,
                edge=True,
                pad=pad,
                pad_top=pad_top,
                nowrap=True,
                **kwargs
            )

        # List of commands
        else:
            return output_text(
                openad_help.queried_commands(matching_commands, inp=inp, starts_with_only=starts_with_only),
                self,
                pad=pad,
                pad_top=pad_top,
                nowrap=True,
                **kwargs
            )

    def preloop(self):
        """CMD class called function: Preloop is called by cmd to get an update the history file each History File"""
        if readline and os.path.exists(self.histfile):
            # note history files can get corrupted so using try to compensate
            try:
                readline.read_history_file(self.histfile)
            except Exception:  # pylint: disable=broad-exception-caught # do not need to know exception
                # Create history file in case it doesn't exist yet.
                # - - -
                # To trigger:
                # >> create new workspace foobar
                # >> ctrl+c
                # (Reboot)
                readline.write_history_file(self.histfile)

    def postloop(self):
        """CMD class called function: Post loop is called by cmd to get an update the history file"""
        readline.set_history_length(self.histfile_size)
        readline.write_history_file(self.histfile)

    def add_history(self, inp):
        """CMD class called function: adds history file"""
        readline.add_history(inp)

    def complete(self, text, state):
        """CMD class called function:
        This is the Auto Complete method that gets called on Forward Tab.
        Currently parsing list pyparsing statements and finding the statement that it fails against at a character
        furthest along the command the statement string is the method used.
        This is an area for improvement further along the line

        Note: this will only match up until optional components in the command,
          Parser does not predict options uptake in command
        """
        if state == 0:
            orig_line = readline.get_line_buffer()
        i_s = 0
        yy = []

        if len(orig_line.split()) > 1:
            orig_word = orig_line.split()[len(orig_line.split()) - 1]
        else:
            orig_word = orig_line

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

                x = x.replace(orig_line, "")
                if (
                    x.split(",")[0].find("Expected CaselessKeyword") > -1
                    or x.split(",")[0].find("Expected Keyword") > -1
                    or x.split(",")[0].find("Expected {") > -1
                ) and x.split(",")[0].find("'" + orig_word.lower()) > -1:
                    if len(orig_line.split(">>")) > 1:
                        started_word = str(orig_word).strip()

                        for match in x.split("'"):
                            if match.upper().startswith(started_word.upper().strip()):
                                # print(match.upper() + "     " + started_word.upper())
                                # print(match[len(started_word) :])

                                readline.insert_text(match[len(started_word.strip()) :])
                                readline.redisplay()
                                return ""

                    yy = x.split(",")[0].split("'")[1]

                    readline.insert_text(yy[len(orig_word) :])
                    readline.insert_text(" ")
                    readline.redisplay()  # readline redisplay needed to push Macos prompt to update
                    return ""  # return Nothing Changed

        # Look for a whole word match

        for i in test_list:
            if error_col_grabber(str(i)) < best_fit:  # If worse match continue
                continue

            if len(i) > 1:
                if (
                    x.split(",")[0].find("Expected CaselessKeyword") > -1
                    or x.split(",")[0].find("Expected Keyword") > -1
                    or x.split(",")[0].find("Expected {") > -1
                ) and x.split(",")[1].find("at char 0") > -1:
                    if (
                        str(str(i[1]).split(",", maxsplit=1)[0].split("Keyword")[1].split("'")[1]).strip().upper()
                        == str(i[0] + x.split(",")[0].split("Keyword")[1].split("'")[1]).strip().upper()
                    ):
                        readline.insert_text(x.split(",")[0].split("Keyword")[1].split("'")[1].strip())
                        readline.insert_text(" ")
                        readline.redisplay()

                        return ""
                    if (
                        str(str(i[1]).split(",", maxsplit=1)[0].split("{")[1].split("'")[1]).strip().upper()
                        == str(i[0] + x.split(",")[0].split("{")[1].split("'")[1]).strip().upper()
                    ):
                        readline.insert_text(x.split(",")[0].split("{")[1].split("'")[1].strip())

                        readline.insert_text(" ")
                        readline.redisplay()

                        return ""
                    if len(orig_word.split(">>")) > 1:
                        started_word = str(orig_word.split(">>")[-1]).strip()

                        for match in x.split("'"):
                            if match.upper().startswith(started_word.upper().strip()):
                                # print(match.upper() + "     " + started_word.upper())
                                # print(match[len(started_word) :])

                                readline.insert_text(match[len(started_word.strip()) :])
                                readline.redisplay()
                                return ""
                    continue
        # if failed previously scan look for successfuly space is next logical character

        for i in test_list:
            if error_col_grabber(str(i)) < best_fit:  # If worse match continue
                continue

            if len(i) > 1:
                c = i[1]
                x = c.explain()
                x = x.replace(orig_line, "")

                if (
                    (
                        x.split(",")[0].find("Expected CaselessKeyword") > -1
                        or x.split(",")[0].find("Expected Keyword") > -1
                        or x.split(",")[0].find("Expected {") > -1
                    )
                    and x.split(",")[0].find("'" + orig_word.lower()) == -1
                    or x.split(",")[0].find("'" + orig_word.split(">>").lower()) == -1
                ):
                    spacing = ""
                    if len(orig_line) == len(i[0]):
                        spacing = " "

                    if len(orig_word.split(">>")) > 1:
                        started_word = str(orig_word.split(">>")[-1]).strip()

                        if orig_word.strip()[-1] == ">":
                            if len(x.split("'")) > 2:
                                readline.insert_text(x.split("'")[1])
                                readline.redisplay()

                            return ""

                        for match in x.split("'"):
                            if match.upper().startswith(started_word.upper().strip()):
                                # print(match.upper() + "     " + started_word.upper())
                                # print(match[len(started_word) :])
                                readline.insert_text(match[len(started_word.strip()) :])
                                readline.redisplay()
                                return ""

                    if error_col_grabber(x) - 1 < len(orig_line):
                        if len(orig_line[error_col_grabber(x) - 1 : len(orig_line)].strip()) > 0:
                            return []
                    readline.insert_text(spacing + x.split(",")[0].split("Keyword")[1].split("'")[1].strip())
                    readline.insert_text(" ")
                    readline.redisplay()
                    return ""  # return nothing changed

                if (
                    str(str(i[1]).split(",", maxsplit=1)[0].split("{")[1].split("'")[1]).strip().upper()
                    == str(i[0] + x.split(",")[0].split("{")[1].split("'")[1]).strip().upper()
                ):
                    readline.insert_text(x.split(",")[0].split("{")[1].split("'")[1].strip())
                    readline.insert_text(" ")
                    readline.redisplay()
                    return ""
        # look for a bracket match

        for i in test_list:
            if error_col_grabber(str(i)) < best_fit:
                continue

            if len(i) > 1:
                c = i[1]
                x = c.explain()

                x = x.replace(orig_line, "")

                if x.split(",")[0].find("Expected '('") > -1 or x.split(",")[0].find("Expected ')'") > -1:
                    if x.split(",")[0].find("Expected '('") > -1:
                        readline.insert_text("(")
                    else:
                        readline.insert_text(")")

                    readline.redisplay()
                    return ""  # return Nothing Changed

                if x.split(",")[0].find("Expected string enclosed in '\"'"):
                    if x.find("Expected W:(A-Za-z, #(),-/-9=A-Z_a-z)"):
                        return ""
                    readline.insert_text("'")
                    readline.redisplay()
                    return ""  # return Nothing Changed

                else:
                    pass

        return ""  # return Nothing Changed

    # Catches the exit command
    def do_exit(self, dummy_inp_do_not_remove):
        """CMD Funcion: called on exit command"""
        write_registry(self.settings, self, True)
        delete_session_registry(self.session_id)
        # exiting the application. Shorthand: x q.
        return True

    # prevents on return on a blank line default receiving the previous input to run
    def emptyline(self):
        """prevents on return on a blank line default receiving the previous input to run"""

    def default(self, line):
        """Default method call on hitting of the return Key, it tries to parse and execute the statements."""
        inp = line  # assigning line to input value

        x = None
        if convert(inp).split()[-1] == "?" and (
            not convert(inp).lower().startswith("tell me") or convert(inp).lower() == "tell me ?"
        ):
            # ... ?
            return self.do_help(inp, display_info=False)

        try:
            try:
                self.settings = load_registry(self)
            except Exception as e:  # pylint: disable=broad-exception-caught # error could be unknown
                # Brutal situation where someone hit clear sessions in another session
                # , shut down abruptly so as not to kill registry file.
                output_error(
                    "Fatal error: the session registry is not avaiable, performing emergency shutdown" + str(e),
                    cmd_pointer=self,
                    return_value=False,
                )
                self.do_exit("exit emergency")

            y = self.current_statement_defs.parseString(convert(inp), parseAll=True)
            x = lang_parse(self, y)
            self.prompt = refresh_prompt(self.settings)
            self.memory.before_command()
        except Exception as err1:  # pylint: disable=broad-exception-caught # error could be unknown
            error_descriptor = None
            error_col = -1
            invalid_command = False
            i_s = 0

            while i_s < len(self.current_statements):
                a, b = self.current_statements[i_s].runTests(
                    convert(inp), printResults=False, fullDump=False, parseAll=True
                )

                for i in b:
                    if len(i) > 1:
                        invalid_command = True

                        c = i[1]
                        try:
                            x = c.explain()
                        # we do not know what the error could be, so no point in being more specific
                        except Exception as err:  # pylint: disable=broad-exception-caught
                            return output_error(msg("err_unknown", err1, split=True), self, return_val=False)

                        if x.find("Expected CaselessKeyword") > -1 and x.find("at char 0") == -1:
                            if error_col < error_col_grabber(x):
                                error_descriptor = x.replace("CaselessKeyword", "keyword").replace(
                                    "ParseException:", "Syntax Error:: "
                                )
                                error_col = error_col_grabber(x)

                        elif x.find("found end of text") > -1 and x.find("at char 0") == -1:
                            if error_col < error_col_grabber(x):
                                error_descriptor = x.replace("ParseException:", "Syntax Error:: ")
                                error_col = error_col_grabber(x)

                        else:
                            if error_col < error_col_grabber(x):
                                error_descriptor = x
                                error_col = error_col_grabber(x)
                i_s = i_s + 1

            # Print error
            if invalid_command:
                # Example input: `search for molecules in parents`

                # To be double checked but... this is an impossible condition. value will always be 1 or more.

                if error_col_grabber(error_descriptor) == 0:
                    return output_error(msg("err_invalid_cmd", strip_tags(msg("run_?")), split=True), self, pad=0)
                else:
                    # Determine if the user input is a partially correct command
                    # or an incorrect command.
                    # Try to match the users intended command and display its help, or
                    # suggest possible commands they were trying to type.

                    # Isolate part of the error message with command & arrow.
                    error_msg = error_descriptor.split("Syntax")[0].splitlines()
                    error_msg = error_msg[0] + "\n" + error_msg[1]

                    # Isolate the string we want to use to search for related commands.
                    if error_col_grabber(error_descriptor) == 1:
                        # No full word recognized --> grab first word.
                        help_ref = error_first_word_grabber(error_descriptor)
                    else:
                        # One of more words recognized --> grab all recognized words.
                        help_ref = inp[0 : error_col_grabber(error_descriptor) - 1]

                    # Not for printing
                    # Fetch commands matching the entire input.
                    # Example input -> `search for molecules in parents`
                    do_help_output_A = self.do_help(
                        inp + " ?", return_val=True, jup_return_format="plain", starts_with_only=True
                    )

                    # Not for printing
                    # Fetch commands matching recognized words, plus
                    # the first letter of the first unrecognized word.
                    # Example input -> `search for molecules in p`
                    do_help_output_B = self.do_help(
                        inp[0 : error_col_grabber(error_descriptor)] + " ?",
                        jup_return_format="plain",
                        starts_with_only=True,
                        return_val=True,
                    )

                    # Not for printing
                    # Fetch commands matching recognized words, or the first word.
                    # Example input -> `search for molecules in`
                    do_help_output_C = self.do_help(
                        help_ref + " ?", return_val=True, jup_return_format="plain", starts_with_only=True
                    )

                    # Check for scenario A, B, C in that order.
                    if "No commands" not in do_help_output_A:
                        # Scenario A
                        show_suggestions = True
                        multiple_suggestions = "Commands starting with" in str(do_help_output_A)
                        help_ref = inp
                    elif "No commands" not in do_help_output_B:
                        # Scenario B
                        show_suggestions = True
                        multiple_suggestions = "Commands starting with" in do_help_output_B
                        help_ref = inp[0 : error_col_grabber(error_descriptor)]
                    else:
                        # Scenario C
                        show_suggestions = "No commands" not in do_help_output_C
                        multiple_suggestions = "Commands starting with" in do_help_output_C

                    # If there are still no suggestions, we loop backwards through
                    # the input string letter by letter until something pops up.
                    if not show_suggestions:
                        error_col = error_col_grabber(error_descriptor)
                        while error_col > 1 and not show_suggestions:
                            error_col = error_col - 1
                            # Not for printing
                            do_help_output_A = self.do_help(
                                inp[0:error_col] + " ?",
                                return_val=True,
                                jup_return_format="plain",
                                starts_with_only=True,
                            )
                            show_suggestions = "No commands" not in str(do_help_output_A)
                            multiple_suggestions = "Commands starting with" in str(do_help_output_A)
                            help_ref = inp[0:error_col]

                    # Display error.
                    output_error(msg("err_invalid_cmd", error_msg, split=True), self, return_val=False)
                    if show_suggestions:
                        if not multiple_suggestions:
                            output_text("<yellow>You may want to try:</yellow>", self, return_val=False)
                            pad_top = 1  # Single command should get a linebreak before and after.
                        else:
                            pad_top = 0  # List of commands should not get padded.

                        # Example to trigger this: `list xxx`
                        self.do_help(help_ref + " ?", starts_with_only=True, return_val=False, pad_top=pad_top)
                        output_text(msg("run_?"), self, return_val=False, pad=1)
                    return

            else:
                output_error(msg("err_unknown"), self, return_val=False)
                return

        if self.refresh_train is True:
            output_train_statements(self)
            self.refresh_train = False
        if self.notebook_mode is True:
            return x
        elif self.api_mode is False:
            if x is not None and not isinstance(x, bool):
                print(x)
            else:
                return


# Returns the error positioning in the statement that has been parsed.
def error_col_grabber(error):
    """Returns the error positioning in the statement that has been parsed."""
    e = error.split("col:")[1]
    e1 = e.replace(")", "")
    return int(e1)


def error_first_word_grabber(error):
    """Retuns the error positioning WORD in the statement that has been parsed."""
    word = error.split("found ")[1].split("'")[1]

    return str(word)


# Main execution application
# If the application is called with parameters, it executes the parameters.
# If called without parameters, the command line enters the shell environment.
# History is only kept for commands executed once in the shell.
def api_remote(
    inp: str,
    api_context: dict = {"workspace": None, "toolkit": None},
    api_var_list={},
):
    """
    Receives Notebook Magic Command API calls
    Note: As we do not have a continuous session in notebooks like on a command line.
    The api_context subsititues for this and the in memory persistent toolkit handle cache
    this makes each command smoother and faster after each Magic command request.
    So api_context handles these:
    - What workspace is current
    - What toolkit is in context
    - Caching of all handles
    - It is deliberate that the whole RUNCMD class object is not kept alive as there is no logical
      exit point for magic commands, unlike a command line.
    """

    initialise()

    arguments = inp.split()
    inp = ""  # reset input after splitting into arguments
    a_space = ""  # reset a_space

    # setup for notebook mode
    magic_prompt = RUNCMD(notebook=True)
    magic_prompt.notebook_mode = True

    if api_context["workspace"] is None:
        api_context["workspace"] = magic_prompt.settings["workspace"]
    else:
        x = {"Workspace_Name": api_context["workspace"]}
        set_workspace(magic_prompt, x)

    if api_context["toolkit"] is None:
        api_context["toolkit"] = magic_prompt.settings["context"]
    else:
        x = {"toolkit_name": api_context["toolkit"]}
        set_context(magic_prompt, x)

    magic_prompt.api_variables = api_var_list
    # We now manage history. The history sometimes gets corrupted through no fault of ours.
    # If so, we just reset it.
    try:
        readline.read_history_file(magic_prompt.histfile)
    except Exception:  # pylint: disable=broad-exception-caught # could be a number of errors
        readline.add_history("")
        readline.write_history_file(magic_prompt.histfile)
        readline.read_history_file(magic_prompt.histfile)
    for i in arguments:
        inp = inp + a_space + i
        a_space = " "

    # Check for Help from a command line request
    if len(inp.strip()) > 0:
        if (
            inp.split()[0] == "?"
            or (
                inp.split()[-1] == "?"
                and (not convert(inp).lower().startswith("tell me") or convert(inp).lower() == "tell me ?")
            )
            or inp.strip() == "??"
        ):
            if inp.strip() == "?":
                inp = ""
            elif inp.strip() == "??":
                inp = "?"

            # return magic_prompt.do_help(inp, display_info=True)
            # Triggered by magic commands, eg. `%openad list files ?`
            # display_info = inp.split()[0] == "?" and inp.strip() != "??"
            # return magic_prompt.do_help(inp.strip(), display_info=display_info)

            # Triggered by magic commands, eg. `%openad ? list files`
            starts_with_qmark = len(inp) > 0 and inp.split()[0] == "?" and inp.strip() != "??"
            return magic_prompt.do_help(inp.strip(), display_info=starts_with_qmark)

        # If there is a argument and it is not a help attempt to run the command.
        # Note, may be possible add code completion here #revisit
        else:
            magic_prompt.preloop()
            magic_prompt.add_history(inp)
            magic_prompt.postloop()
            readline.write_history_file(magic_prompt.histfile)

            result = magic_prompt.default(inp)
            api_context["workspace"] = magic_prompt.settings["workspace"]
            api_context["toolkit"] = magic_prompt.settings["context"]

            magic_prompt.do_exit("dummy do not remove")
            if result is not True and result is not False:
                return result


def cmd_line():
    """Entry point for command line interface"""
    initialise()
    inp = ""
    a_space = ""

    for i in sys.argv[1:]:
        inp = inp + a_space + i
        a_space = " "
    try:
        command_line = RUNCMD()
    except KeyboardInterrupt:
        output_error("Keyboard Initiated Exit before OpenAD Initialised")
        return

    # Check for help from a command line request.
    if len(inp.strip()) > 0:
        words = inp.split()
        if inp.split()[0] == "?" or inp.split()[-1] == "?" or inp.strip() == "??":
            # Impossible clause...?
            # you can't add `?` to `openad <command>` because the terminal interprets the `?` separately.
            if inp.strip() == "?":
                inp = ""
            elif inp.strip() == "??":
                inp = "?"

            # Triggered by running commands from main terminal prepended with `openad`.
            starts_with_qmark = len(inp) > 0 and inp.split()[0] == "?" and inp.strip() != "??"
            command_line.do_help(inp.strip(), display_info=starts_with_qmark)

        # If user wants to run command line and specify toolkit, for a specific command:
        elif words[0].lower() == "-s" and len(words) > 3:
            set_workspace(command_line, {"Workspace_Name": words[1].upper()})
            set_context(command_line, {"toolkit_name": words[2].upper()})
            if (
                command_line.settings["workspace"] == words[1].upper()
                and command_line.settings["context"] == words[2].upper()
            ):
                command_line.preloop()
                command_line.add_history(str(" ".join(words[3:])).strip())
                command_line.postloop()
                command_line.default(str(" ".join(words[3:])).strip())
        else:
            # If there is an argument and it is not help, attempt to run the command
            # Note, may be possible add code completion here #revisit
            command_line.preloop()
            command_line.add_history(inp.strip())
            command_line.postloop()
            command_line.default(inp.strip())
        command_line.do_exit("dummy do not remove")
    else:
        # If no argument passed then enter
        lets_exit = False
        while lets_exit is False:
            try:
                # The cmdloop parameter controls the startup screen, it overrides self.intro.
                command_line.cmdloop(splash(command_line.settings["context"], command_line, startup=True))
                lets_exit = True
            except KeyboardInterrupt:
                command_line.postloop()
                if confirm_prompt("Are you sure you wish to exit?", default=True):
                    lets_exit = True
                    command_line.do_exit("dummy do not remove")
            except Exception as err:  # pylint: disable=broad-exception-caught
                # we do not know what the error could be, so no point in being more specific
                output_error(msg("err_invalid_cmd", err, split=True), command_line)


if __name__ == "__main__":
    cmd_line()
