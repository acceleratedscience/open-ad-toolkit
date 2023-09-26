from ad4e_opentoolkit.helpers.general import singular, is_toolkit_installed
from ad4e_opentoolkit.helpers.output import msg, output_text, output_error
# Importing our own plugins.
# This is temporary until every plugin is available as a public pypi package.
from ad4e_opentoolkit.plugins.style_parser import style

line_break = "<br>"
heading = "<h3>"
heading_end = "</h3>"
bold_begin = "<b>"
bold_end = "</b>"
italic_begin = "<i>"
italic_end = "</i>"
spacing = " "
six_spacing = "      "
bullet = " - "


# Create the help dictionary object for a command.
def help_dict_create(name: str, command: str, description: str, url: str = None, category: str = 'Uncategorized'):
    return {
        'category': category,
        'name': name,
        'command': command,
        'description': description,
        'url': url,
    }


def all_commands(commands: list, toolkit_name: str = None, toolkit_current: object = None, cmd_pointer: object = None):
    """
    Return xml string listing all available commands organized by category.

    Command: `?`
    """
    def _compile(commands, toolkit_name=None):
        commands_organized = {}

        # Organize commands by category.
        for i in commands:
            command = i['command']
            try:
                category = i['category']
            except BaseException:
                category = 'Uncategorized'
            if category in commands_organized:
                commands_organized[category].append(command)
            else:
                commands_organized[category] = [command]

        # Compile text.
        output = [f'<h1>Available Commands - {toolkit_name if toolkit_name else "Main"}</h1>']
        if toolkit_name and not is_toolkit_installed(toolkit_name, cmd_pointer):
            err_msg = output_error(
                msg('fail_toolkit_not_installed', toolkit_name, split=True),
                cmd_pointer, return_val=True, nowrap=True)
            output.append(err_msg)
        elif len(commands_organized):
            output.append('')

        if len(commands_organized):
            for category, commands in commands_organized.items():
                output.append(f"{category}:")
                for command in commands:
                    output.append(f'<cmd>{command}</cmd>')
                output.append('')
        else:
            output.append('<error>No commands found.</error>')
        return '\n'.join(output)

    #
    #

    if toolkit_name:
        return _compile(commands, toolkit_name)
    else:
        main_commands = _compile(commands)
        if toolkit_current:
            toolkit_commands = '\n\n\n' + _compile(toolkit_current.methods_help, toolkit_current.toolkit_name)
        else:
            toolkit_commands = ''
        return main_commands + toolkit_commands


def queried_commands(command_results: list, inp: str = None, query_type: str = None):
    """
    Return xml string listing all commands that match the query.


    Command: `list ?` --> 'List' matches a command word, so we
                          display all commands using this word.
    Command: `lis ?`  --> 'Lis' doesn't match a command word, so
                          we display all commands startig with "lis".
    """
    import re

    # Compile title.
    if query_type == 'word_match':
        inp = singular(inp)
        title = f'Commands containing "{inp}"'
    elif query_type == 'starts_with':
        title = f'Commands starting with "{inp}..."'

    # Compile text.
    output = [f'<yellow>{title}</yellow>']
    for command in command_results:
        command_str = command['command']
        if query_type == 'word_match':
            command_str = re.sub(fr'(?<!<){inp}(s?)(?![^<>]*?>)', fr'<underline>{inp}\1</underline>', command_str)
        output.append(f"- <cmd>{command_str}</cmd>")
    return '\n'.join(output)


def command_details(command: list):
    """
    Return xml string listing listing a single command with its description.

    Command: `<command> ?`
    """

    paragraph_width = 100

    # Style command
    the_command = style(f"<cmd>{command['command']}</cmd>", width=paragraph_width)

    # Separator
    sep_len = min(len(command['command']), paragraph_width)
    sep = '<soft>' + sep_len * '-' + '</soft>'

    # Style description
    description = style(command['description'], width=paragraph_width)

    return '\n'.join([the_command, sep, description])


# Display advanced help
def advanced_help():
    return '<warning>Advanced help is yet to be implemented.</warning>'


class openad_help():
    help_orig = []
    help_current = []

    def add_help(self, help_dict):
        for i in help_dict:
            self.help_current.append[i]

    def reset_help(self):
        self.help_current = self.help_orig.copy()
