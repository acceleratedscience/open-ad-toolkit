import sys
import os
import re
import readline
from IPython.display import display

from helpers.output import msg, output_text, output_error
from helpers.style_parser import strip_tags


# Refreshes the command prompt when in the shell.
def refresh_prompt(settings):
    if settings['context'] is not None:
        #prompt = ' \u001b[7m ' + settings['context'] + ' \u001b[0m '  # Reverse & reset
        prompt = settings['context']+'->'   # Reverse & reset
    else:
        prompt = 'ADCCL:'
    if settings['workspace'] is not None:
        prompt = prompt + settings['workspace']
    prompt = prompt + ' >>  '
    return prompt


# Todo: also check for API
#


def is_notebook_mode():
    """ Return True if we are running inside a Jupyter Notebook or Jupyter Lab. """
    try:
        get_ipython()
        return True
    except BaseException:
        return False


def remove_lines(count=1):
    """ Remove the last printed line(s) from the CLI. """
    if is_notebook_mode():
        # Jupyter
        # In Jupyter you can't clear a single line, only the entire cell output.
        from IPython.display import clear_output
        clear_output(wait=True)
    else:
        # CLI
        while count > 0:
            count -= 1
            sys.stdout.write('\033[F')  # Move the cursor up one line
            sys.stdout.write('\033[K')  # Clear the line
            sys.stdout.flush()  # Flush the output buffer


def convertTuple(tup):
    # initialize an empty string
    if isinstance(tup, tuple):
        str = ''
        space = ''
        for item in tup:
            str = str + item
        return str
    else:
        return tup


def singular(string):
    return re.sub(r's$', '', string)


def parse_path_tree(path_string):
    # Normalize the path string to use the appropriate separator for the current system.
    path = os.path.normpath(path_string)

    # Cut off first slash if it exists.
    if (path[0] == '/'):
        path = path[1:]

    # Remove file name and return the tree
    tree = path.split('/')[:-1]
    return tree


# Confirm promt for True or False Questions
def confirm_prompt(question: str) -> bool:
    reply = None
    while reply not in ('y', 'n'):
        try:
            question_formatted = output_text(f'<yellow>\n{question}</yellow>', return_val=True, pad=1)
            if is_notebook_mode():
                display(question_formatted)
                reply = input('(y/n): ').casefold()
            else:
                space = '' if question_formatted[-1] == '\n' else ' '
                reply = input(question_formatted + f'{space}(y/n): ').casefold()
            readline.remove_history_item(readline.get_current_history_length() - 1)
        except KeyboardInterrupt:
            print('')
            return True
    if reply == 'y':
        print('')
        return True


# Return boolean and formatted error message if other sessions exist.
def other_sessions_exist(cmd_pointer):
    from global_var_lib import _meta_registry_session as _meta_registry_session
    file_list = os.listdir(os.path.dirname(_meta_registry_session))
    try:
        file_list.remove('registry.pkl' + cmd_pointer.session_id)
    except BaseException:
        pass

    if len(file_list) > 0:
        return True, output_error(msg('abort_clear_sessions', split=True), cmd_pointer, return_val=False)
    else:
        return False, None


# Return user input.
def user_input(cmd_pointer, question):
    """
    Basically the same as input(), but with some extra styling and history disabled.
    """
    prompt = output_text(f'<yellow>{question}: </yellow>', cmd_pointer, return_val=True, jup_return_format='plain')
    text = input(prompt)
    return text


# Return list of available toolkit names.
def get_toolkits():
    import os
    
 
    folder_path = os.path.dirname(os.path.abspath(__file__))+'/../user_toolkits'  # Replace 'foo' with the actual path to the folder
    toolkit_names = [name.upper() for name in os.listdir(folder_path) if os.path.isdir(os.path.join(folder_path, name))]
    return toolkit_names


# Return boolean if toolkit is installed.
def is_toolkit_installed(toolkit_name, cmd_pointer=None):
    return cmd_pointer and toolkit_name and toolkit_name.upper() in cmd_pointer.settings['toolkits']

#
#
#
#
#
#
#
#
# NOT USED
#
#
#
#
#
#
#
#

# NOT USED
# Operating system


def get_platform():
    platforms = {
        'linux1': 'Linux',
        'linux2': 'Linux',
        'darwin': 'OS X',
        'win32': 'Windows'
    }
    if sys.platform not in platforms:
        return sys.platform
    return platforms[sys.platform]


# NOT USED
# Alternative to spinner that exits more gracefully on ctrl+c.
# Abandoning this for now because it's blocking the process.
def loader(text='', anim=["◐", "◓", "◑", "◒"], no_format=False, exit_msg=None, on_abort=None):
    import asyncio
    interval = 0.1
    sys.stdout.write('\033[?25l')  # Hide cursor

    async def _loop(text='', i=0, line_length=0):
        # print('>', line_length)
        str = anim[i]
        sys.stdout.write('\r' + ' ' * line_length + '\r')
        sys.stdout.flush()
        # Make text soft gray.
        text_formatted = text if no_format or not text else f'\u001b[90m{text}...\u001b[0m'
        line = '\r' + str + ' ' + text_formatted
        sys.stdout.write(line)  # This should return line_length but it doesn't for some reason.
        line_length = len(line)
        i = (i + 1) % len(anim)
        # time.sleep(interval) # Sync
        # _loop(text, i, line_length) # Sync
        await asyncio.sleep(interval)
        # await asyncio.run(_loop(text, i, line_length))

        try:
            await _loop(text, i, line_length)
            # _loop(text, i, line_length)
        except KeyboardInterrupt:
            print('\b \b \b' + '\b \b' * line_length + '\r')  # Clear the line.
            sys.stdout.write('\033[F')  # Move the cursor up one line
            if exit_msg:
                print(exit_msg)  # Move the cursor up one line
            if on_abort:
                on_abort()

    # _loop(text=text) # Sync
    asyncio.run(_loop(text=text))
    sys.stdout.write('\033[?25h')  # Show cursor


# NOT USED
# Clear the current line.
# Couldn't get this to work, because I can't clear the buffer.
# As a result, (eg.) when you try ctrl-c twice with n answer,
# the second time it will have the first time's response in the buffer,
# causing it to mess up the layout.
def clear_current_line():
    buffer = readline.get_line_buffer()
    line_length = len(buffer)
    readline.clear_history()
    eraser = '\b \b' * line_length
    sys.stdout.write(eraser)
    # readline.insert_text(' ')
    # print(len(buffer), buffer)
