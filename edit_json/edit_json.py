import ast
import os
import re
from functools import reduce
import operator
import json
from blessed import Terminal


class EditJson:
    """
    Edit a JSON file with the CLI.

    Parameters:
    path:   Path to a JSON file. If no path is provided, we will ask for one.
            If no file is found, a new one will be created based on the schema.
    schema: Optional data schema that contains data types and help text per field.
            This throws an error when invalid values are submitted.

    Supported data types:
    bool, str, int (min, max), float (min, max), enum (options, case_sensitive)

    Lists of data types are also supported:
    ['bool'], ['str'], ['int'], ['float'], ['enum']

    All values are optional unless required is set to True.


    ----------------------------------------

    Usage:
    import edit_json

    # Let the user select a JSON file.
    edit_json()

    # Load a JSON file without a schema,
    # and thus no validation or help strings.
    edit_json('settings.json')

    # Load a JSON file and a schema.
    # This will validate the JSON file against
    # the schema and display help strings per field.
    # When no JSON file is found, one will be created
    # based on the schema.
    edit_json('settings.json', 'settings-schema.json')

    # Run a demo with dummy data.
    edit_json(demo=True)

    # Consult the docstring for usage info.
    help(edit_json)

    ----------------------------------------

    Example schema:
    {
        "foo": {
            "type": "bool",
            "required": true,  # Make the field mandatory.
            "default": true,  # When no JSON file is found, a new one will be created with this value.
            "help": "This is a required boolean value."
        },
        "bar": {
            "type": "str",
            "help": "This is a string."
        },
        "baz": {
            "type": "int",
            "min": 0,
            "max": 10,
            "help": "This is an integer between 0 and 10."
        },
        "qux": {
            "type": "float",
            "min": 0,
            "help": "This is a float above zero."
        }
        "quux": {
            "type": "enum",
            "values": ["foo", "bar", "baz"],
            "case_sensitive": false, # Default is true
            "help": "This is an enum."
        },
        "quuz": {
            "type": ["int"],
            "help": "This is a list of integers."
        },
    }
    """
    term = Terminal()
    colors = {
        'white': (255, 255, 255),
        'soft_white': (150, 150, 150),
        'black': (0, 0, 0),
        'soft_black': (90, 90, 90),
        'blue': (50, 230, 240),
        'blue_soft': (21, 69, 71),
        'yellow': (255, 255, 0),
        'red': (255, 0, 0)
    }
    styles = {
        'reg': term.on_color_rgb(*colors['black']) + term.color_rgb(*colors['white']),
        'info': term.on_color_rgb(*colors['black']) + term.cyan,
        'info_highlight': term.on_color_rgb(*colors['blue_soft']) + term.cyan,
        'info_inverse': term.on_cyan + term.color_rgb(*colors['black']),
        'soft': term.on_color_rgb(*colors['black']) + term.color_rgb(*colors['soft_white']),
        'sel': term.on_color_rgb(*colors['white']) + term.color_rgb(*colors['black']),
        'edit': term.on_color_rgb(*colors['soft_black']) + term.color_rgb(*colors['white']),
        'cursor': term.on_color_rgb(*colors['white']) + term.color_rgb(*colors['black']),
        'indent': term.on_color_rgb(*colors['black']) + term.color_rgb(*colors['soft_black']),
        'help': term.on_color_rgb(*colors['black']) + term.color_rgb(*colors['yellow']),
        'error': term.on_red + term.color_rgb(*colors['white']),
    }

    def __init__(self):
        # Initialize
        self.log_value = ''  # For debugging using self.log('foo')
        self.json_file_path = ''
        self.json_file_name = ''
        self.data = {}  # To store JSON data until it gets saved to file.
        self.edit_mode = False
        self.text_input = ''  # To store typing until you press enter.
        self.cursor_pos = 0
        self.input_before_cursor = ''
        self.input_after_cursor = ''
        self.error = ''  # General error message
        self.inline_errors = {}  # Error messages per field
        self.illegal_errors = {}  # When a field is not present in the schema.
        self.selected_index = 0
        self.selected_key = ''
        self.selected_value = ''
        self.selected_help = ''
        self.listen_for_action = False  # When you type : we wait for 'w' or 'q' command
        self.left_col_width = 20
        self.show_help = False
        self.show_confirm_exit = False
        self.indent = '|   '  # Indentation per level
        self.paths = []  # See store_path()
        self.path = None
        self.schema = None
        # Lets us hide the file path of demo.json & block file saving.
        self.demo = False

        self.debug = True  # Display logger

    def __call__(self, path=None, schema=None, demo=False):
        self.__init__()

        # Load demo data if requested.
        if demo:
            dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'demo')  # noqa
            path = os.path.join(dir, 'demo.json')
            schema = os.path.join(dir, 'demo-schema.json')
            self.demo = True

        # Load schema data.
        if schema:
            try:
                with open(schema, encoding='UTF-8') as file:
                    # raise FileNotFoundError("Test.")  # To simulate error.
                    self.schema = json.load(file)
            except (FileNotFoundError):
                print(self.styles['error'],
                      f'Invalid schema file path: {schema}', self.styles['reg'])
            except (json.JSONDecodeError):
                print(self.styles['error'],
                      f'Schema file is not a valid JSON: {schema}', self.styles['reg'])

        # Abort if there's a problem with the schema file.
        if schema and not self.schema:
            return

        # Load JSON data.
        abort = False
        while True:
            self.json_file_path = path if path else input(
                'Enter path to JSON file: ')
            try:
                with open(self.json_file_path, encoding='UTF-8') as file:
                    self.data = json.load(file)
                    self.json_file_name = file.name
                    break
            except (FileNotFoundError):
                if not path:
                    print(self.styles['error'],
                          f'File not found: {self.json_file_path}', self.styles['reg'])
                    abort = True
                else:
                    # When no JSON file is found at the provided path, we create one.
                    if not self.data and schema:
                        self.create_file_from_schema()
                    else:
                        print(self.styles['error'],
                              f'File not found: {self.json_file_path}', self.styles['reg'])
                        abort = True
            except (json.JSONDecodeError):
                print(self.styles['error'],
                      f'Invalid JSON file: {self.json_file_path}', self.styles['reg'])
                abort = True
            finally:
                # Break is needed for when we're asking
                # the user for a path using the while loop.
                if path:
                    break

        # Abort if there's a problem with the JSON file.
        if abort:
            return

        # Parse JSON data to set up UI.
        self.store_paths()
        self.set_left_col_width()

        # Validate the data & render the UI.
        self.validate_all_values()
        self.render_list()

    # For debugging: print debug data above the JSON content.
    def log(self, v):
        self.log_value = v

    # Create JSON file from schema.
    def create_file_from_schema(self):

        def add_one_level(schema, data):
            for key, value in schema.items():
                if type(value) == dict and not ('type' in value or 'help' in value):
                    data[key] = {}
                    add_one_level(value, data[key])
                else:
                    data[key] = value['default'] if 'default' in value else ''

        self.data = {}
        add_one_level(self.schema, self.data)

    # Store path
    def store_paths(self):
        """
        Store a flat array with all nested keys as tuples.

        This is used to navigate between fields when pressing UP/DOWN keys.

        Data:
        {
            foo: 1
            bar: {
                a: 2,
                b: 3
            }
        }

        Paths:
        [('foo'), ('bar', 'a'), ('bar', 'b')]
        """
        self.paths = []

        def store_one_level(data, parent_keys=tuple()):
            for i, (key, value) in enumerate(data.items()):
                is_dict = isinstance(value, dict)
                if is_dict:
                    store_one_level(value, parent_keys=parent_keys + (key,))
                else:
                    self.paths.append(parent_keys + (key,))

        store_one_level(self.data)

    # Measure the length of all keys on all levels,
    # accounting for indentation, and set the left
    # column width so it fits everything, but cap
    # it at 50 characters.
    def set_left_col_width(self):
        def measure_keys_one_level(data, level=0):
            max_key_length = 0
            for i, (key, value) in enumerate(data.items()):
                is_dict = isinstance(value, dict)
                if is_dict:
                    max_key_length = max(
                        measure_keys_one_level(value, level=level+1),
                        len(key + self.indent * level),
                        max_key_length)
                else:
                    max_key_length = max(
                        len(key + self.indent * level),
                        max_key_length)
            return max_key_length

        max_key_length = measure_keys_one_level(self.data)
        self.left_col_width = min(max_key_length + 2, 50)

    # Render interactive list of key-value pairs.
    def render_list(self):
        """
        Render interactive list of key-value pairs.
        """

        # Display instructions on top.
        def print_instructions():
            print(self.styles['info_inverse'] + ' ▲ ' + self.styles['info'] + ' ' + self.styles['info_inverse'] + ' ▼ ' + self.styles['info'] + ' Navigate   ', end='')  # noqa
            print(self.styles['info_inverse'] + ' ▶ ' + self.styles['info'] + ' Property info   ', end='')  # noqa
            print(self.styles['info_inverse'] + ' :w ' + self.styles['info'] + ' Save   ', end='')  # noqa
            print(self.styles['info_inverse'] + ' :q ' + self.styles['info'] + ' Exit' + self.styles['reg'])  # noqa

        exit_msg = None
        with self.term.fullscreen(), self.term.cbreak(), self.term.hidden_cursor():
            while True:
                # Clear the console.
                print(self.term.home + self.term.clear)

                # Store current value.
                self.path = self.paths[self.selected_index]
                self.selected_key = self.path[-1]
                self.selected_value = self.value_from_path()
                if self.schema:
                    schema_prop = self.value_from_path(
                        data=self.schema)
                    help_str = schema_prop['help'] if schema_prop and 'help' in schema_prop else ''
                    type_str = schema_prop['type'] if schema_prop and 'type' in schema_prop else ''
                    range_str = ''
                    if type_str == 'int' or type_str == 'float':
                        if 'min' in schema_prop and 'max' in schema_prop:
                            range_str = f" {schema_prop['min']}-{schema_prop['max']}"
                        elif 'min' in schema_prop:
                            range_str = f" >{schema_prop['min']}"
                        elif 'max' in schema_prop:
                            range_str = f" <{schema_prop['max']}"
                    self.selected_help = f'({type_str}{range_str}) {help_str}' if type_str else help_str

                # Debugging display.
                # print('Edit mode:', self.edit_mode)
                # print('Cursor:', self.cursor_pos)
                # print('Listen:', self.listen_for_action)
                # print('Input:', self.text_input)
                # print(
                #     'Input:', f'{self.input_before_cursor}/{self.input_after_cursor}')
                # print('Index:', self.selected_index)
                # print('Path:', self.path)
                # print('Selected key:', self.selected_key)
                # print('Selected value:', self.selected_value)
                # print('Selected help:', self.selected_help)
                # print('Paths:', '\n  ' + '\n  '.join(['{}: {}'.format(i, v)
                #       for i, v in enumerate(self.paths)]))
                # print('Inline errors: ', end='')
                # pprint(self.inline_errors)
                # print('\n')

                # Print instruction header.
                print(self.styles['info'] + 'You are editing ' +
                      self.styles['info_highlight'] + f" {'demo.json' if self.demo else self.json_file_name} " + self.styles['reg'])
                if self.listen_for_action:
                    print(':')
                elif self.error:
                    print(self.styles['error'] +
                          self.error + self.styles['reg'])
                elif self.illegal_errors:
                    error_msg = 'The selected JSON file structure does not match the provided schema file. Type :f to remove all illegal values.'
                    print(self.styles['error'] +
                          error_msg + self.styles['reg'])
                else:
                    print_instructions()

                print(self.styles['reg'] + '------')

                if self.show_confirm_exit:
                    print(
                        self.styles['reg'] + 'Are you sure you want to discard your changes?\ny/n')
                else:
                    # Print list.
                    self.print_list()
                    if (self.debug):
                        print('\nLog: --', self.log_value)  # Debug

                # Process user input (ks = keystroke)
                # Note: the while loop waits for the inkey input to be released.
                ks = self.term.inkey(timeout=3, esc_delay=0)
                exit_msg = self.process_key_stroke(ks)
                if exit_msg:
                    break

        if exit_msg != True:
            print(exit_msg)

    # Print list with each keyboard stroke.
    def print_list(self):
        index = 0

        def print_one_level(data, index, indent_level=0, path=()):
            for i, (key, value) in enumerate(data.items()):
                is_dict = isinstance(value, dict)
                indent = self.styles['indent'] + (self.indent * indent_level)
                current_path = path + (key,)
                if is_dict:
                    # Expand one level deeper.
                    print(indent + self.styles['soft'] + key + ':')
                    index = print_one_level(
                        value, index, indent_level=indent_level+1, path=current_path)
                else:
                    is_selected = self.path == current_path
                    illegal_flag = self.styles['error'] + ' ILLEGAL ' + \
                        self.styles['reg'] + \
                        ' ' if current_path in self.illegal_errors else ''

                    # # Display index for debugging:
                    # if not illegal_flag:
                    #     print(f'{index}:', end='')

                    if is_selected and not illegal_flag:
                        # Selected row.
                        print(indent + self.styles['reg'] +
                              f"{(key + ':'):<{self.left_col_width - len(self.indent * indent_level)}}", end='')
                        if self.edit_mode:
                            # Edit mode.
                            input_after_cursor = f'{self.input_after_cursor} '
                            print(self.styles['edit'] + ' ' + self.input_before_cursor +
                                  self.styles['cursor'] + input_after_cursor[:1] +
                                  self.styles['edit'] + input_after_cursor[1:] + self.styles['soft'])
                        else:
                            print(self.styles['sel'] +
                                  f' {value} ' + self.styles['help'] + (' ▶' if self.selected_help and not self.show_help else '') + self.styles['reg'])

                        # Helper text.
                        if self.show_help and self.selected_help:
                            help_str = self.selected_help.replace(
                                '\n', '\n' + (' ' * self.left_col_width) + ' ')
                            print(self.styles['soft'] + (' ' * self.left_col_width) + ' ' +
                                  self.styles['help'] + help_str + self.styles['reg'])
                    else:
                        # All other rows.
                        print(
                            indent + self.styles['soft'] + f"{key + ':':<{self.left_col_width - len(self.indent * indent_level)}} {illegal_flag + self.styles['reg'] + str(value)}")

                    if not illegal_flag:
                        # Inline error message.
                        if self.inline_errors and self.paths[index] in self.inline_errors:
                            print(self.styles['soft'] + (' ' * self.left_col_width) +
                                  self.styles['error'] + f" {self.inline_errors[self.paths[index]]} " + self.styles['soft'])

                        index += 1

            return index

        print_one_level(self.data, index)

    # Recursively loop through data structure and select value based on path.
    def value_from_path(self, data=None, path=None):
        data = data if data else self.data
        path = path if path else self.path
        return reduce(lambda d, k: d.get(k) if isinstance(d, dict) else d[int(k),], path, data)

    # Recursively loop through data structure and set value based on path.
    def value_to_path(self, value, data=None, path=None):
        data = data if data else self.data
        path = path if path else self.path
        reduce(operator.getitem, path[:-1], data)[path[-1]] = value

    # Recursively loop through data structure and remove key based on path.
    def remove_key_from_path(self, data=None, path=None):
        data = data if data else self.data
        path = path if path else self.path
        reduce(operator.getitem, path[:-1], data).pop(path[-1])

    # Keyboard functionality.
    def process_key_stroke(self, ks):
        self.error = ''

        # Confirm screens
        if self.show_confirm_exit:
            if ks == 'y' or ks.name == 'KEY_ENTER':
                self.show_confirm_exit = False
                return 'Your changes havs been discarded.'
            elif ks == 'n' or ks.name == 'KEY_ESCAPE':
                self.show_confirm_exit = False
                return
            else:
                return

        # Listen for actions, or deactivate listener.
        if self.listen_for_action:
            self.listen_for_action = False
            actions = {
                'q': self.confirm_exit,
                'w': self.update_file,
                'f': self.format_file,
            }
            if ks in actions:
                exit_msg = actions[ks]()
                return exit_msg

        if self.edit_mode:
            # Edit mode

            # Deactivate action listener.
            self.listen_for_action = False

            # Map key
            # Keys agnostic of edit mode.
            if ks.name == 'KEY_UP':
                self.store_value()
                self.go_to_prev_field()
            elif ks.name == 'KEY_DOWN':
                self.store_value()
                self.go_to_next_field()
            elif ks.name == 'KEY_ENTER':
                self.store_value()
            elif ks.name == 'KEY_ESCAPE':
                self.edit_mode_exit()
            elif ks.name == 'KEY_BACKSPACE':
                self.input_before_cursor = self.text_input[:self.cursor_pos-1]
                self.input_after_cursor = self.text_input[self.cursor_pos:]
                self.text_input = self.input_before_cursor + self.input_after_cursor
                self.cursor_pos -= 1
            elif ks.name == 'KEY_DELETE':
                self.input_before_cursor = self.text_input[:self.cursor_pos]
                self.input_after_cursor = self.text_input[self.cursor_pos:][1:]
                self.text_input = self.input_before_cursor + self.input_after_cursor
            elif ks.name == 'KEY_LEFT':
                if self.cursor_pos > 0:
                    self.set_cursor(-1)
            elif ks.name == 'KEY_RIGHT':
                if self.cursor_pos < len(self.text_input):
                    self.set_cursor(1)
            elif ks and not ks.is_sequence:
                # Type
                self.type(ks)
            else:
                pass

        else:
            # Browse mode
            if ks.name == 'KEY_UP':
                self.go_to_prev_field()
            elif ks.name == 'KEY_DOWN':
                self.go_to_next_field()
            elif ks.name == 'KEY_RIGHT':
                self.toggle_help(True)
            elif ks.name == 'KEY_LEFT':
                self.toggle_help(False)
            elif ks.name == 'KEY_ENTER':
                self.edit_mode_enter()
            elif ks.name == 'KEY_ESCAPE':
                self.confirm_exit()
            elif ks.name == 'KEY_BACKSPACE':
                self.edit_mode_enter(clear=True)
            elif ks == ':':
                self.listen_for_action = True
            elif ks and not ks.is_sequence:
                # Type
                self.edit_mode_enter(clear=True)
                self.type(ks)

    # Jump to previous field.
    def go_to_prev_field(self):
        self.selected_index = max(0, self.selected_index - 1)
        self.toggle_help(False)
        self.text_input = ''

    # Jump to next field.
    def go_to_next_field(self):
        self.selected_index = min(
            len(self.paths) - 1, self.selected_index + 1)
        self.toggle_help(False)
        self.text_input = ''

    # Show/hide helper text (if present in schema).
    def toggle_help(self, show=True):
        self.show_help = show

    # Enter edit mode.
    def edit_mode_enter(self, clear=False):
        self.edit_mode = True

        # Remove possible field errors.
        if self.path in self.inline_errors:
            self.inline_errors.pop(self.path)

        if not clear:
            # Transfer text input.
            self.text_input = str(self.selected_value)

            # Set cursor at end of line.
            self.set_cursor()

    # Exit edit mode.
    def edit_mode_exit(self):
        self.validate_value()
        self.edit_mode = False
        self.text_input = ''
        self.input_before_cursor = ''
        self.input_after_cursor = ''

    # Update cursor.
    # Direction can be 1, -1, or when set to none,
    # the cursor is set to the end of the line.
    def set_cursor(self, direction=None):
        if direction:
            self.cursor_pos += direction
        else:
            # Set cursor at end of line
            self.cursor_pos = len(self.text_input)

        self.input_before_cursor = self.text_input[:self.cursor_pos]
        self.input_after_cursor = self.text_input[self.cursor_pos:]

    # Enter one character.
    def type(self, ks):
        self.input_before_cursor = self.text_input[:self.cursor_pos] + ks
        self.input_after_cursor = self.text_input[self.cursor_pos:]
        self.text_input = self.input_before_cursor + self.input_after_cursor
        self.cursor_pos += 1

    # Store the text input in self.data.
    def store_value(self):
        try:
            data_type = self.value_from_path(
                data=self.schema, path=self.path + ('type',))
        except:
            data_type=['str']
        
        # Normalize casing for booleans.
        if data_type == ['bool']:
            self.text_input = re.sub(
                r'(?i)(true|false)', lambda x: x.group().capitalize(), self.text_input)

        # Format lists: mixed quotes, irregular spacing, etc.
        try:
            self.text_input = ast.literal_eval(self.text_input)
        except:
            pass

        # Format booleans (this makes boolean input case insensitive).
        if data_type == 'bool' and type(self.text_input) == str:
            if self.text_input.lower() == 'true':
                self.text_input = True
            elif self.text_input.lower() == 'false':
                self.text_input = False

        # Store value in self.data
        self.value_to_path(self.text_input)
        self.edit_mode_exit()

    # Check if there has been any changes.
    def has_changes(self):
        try:
            with open(self.json_file_path, 'r', encoding='UTF-8') as f:
                file_data = json.load(f)
            if self.data == file_data:
                return False
            else:
                return True
        except:
            return False

    # Confirm discarding of changes & exit.
    def confirm_exit(self):
        if (self.has_changes()):
            self.show_confirm_exit = True
        else:
            return True

    # Cycle through self.data and run validate_value on each value.
    def validate_all_values(self):
        for path in self.paths:
            self.validate_value(path)

        # Remove illegal paths from self.paths.
        # This way they are skipped when cycling through the fields.
        for path in self.illegal_errors:
            path_index = self.paths.index(path)
            del self.paths[path_index]

    # Validate the freshly stored value against the data schema.
    def validate_value(self, path=None):

        # Validate enums.
        def _validate_enum(self, values):
            allowed_values = self.value_from_path(
                data=self.schema, path=path + ('options',))

            # Check if the fields is case sensitive (default true).
            case_sensitive = self.value_from_path(
                data=self.schema, path=path + ('case_sensitive',))
            case_sensitive = bool(case_sensitive)

            # Convert allowed values to lowercase if case insensitive.
            allowed_values_lowercase = [str(value).lower()
                                        for value in allowed_values]
            allowed_values = allowed_values if case_sensitive else allowed_values_lowercase

            for value in values:
                value = str(value)
                value = value if case_sensitive else value.lower()
                if value not in allowed_values:
                    if case_sensitive and value.lower() in allowed_values_lowercase:
                        # Value has inproper case.
                        error_msg = "This field is case sensitive."
                        self.inline_errors[path] = error_msg
                    elif value != '':
                        # Value not in enum.
                        error_msg = f"Expected one of these values: {', '.join(allowed_values)}"
                        self.inline_errors[path] = error_msg

        # Validate numbers.
        def _validate_number(self, values, data_type):
            min = self.value_from_path(
                data=self.schema, path=path + ('min',))
            max = self.value_from_path(
                data=self.schema, path=path + ('max',))
            data_type_name = 'an integer' if data_type == 'int' else 'a floating point number'
            for value in values:
                try:
                    if value == '' or value == 0:
                        # Zero or empty is always valid.
                        valid = True
                    elif data_type == 'float':
                        valid = type(value) == float or type(value) == int
                    elif data_type == 'int':
                        valid = type(value) == int
                except:
                    valid = False

                if not valid:
                    # Not an integer/float.
                    error_msg = f"Expected {data_type_name}."
                    self.inline_errors[path] = error_msg
                elif value != '' and (min or min == 0) and float(value) < min:
                    # Too low.
                    error_msg = f"Value must be greater than {min}."
                    self.inline_errors[path] = error_msg
                elif value != '' and (max or max == 0) and float(value) > max:
                    # Too high.
                    error_msg = f"Value must be less than {max}."
                    self.inline_errors[path] = error_msg

        # Validate booleans.
        def _validate_bool(self, values):
            for value in values:
                if value != '' and type(value) != bool:
                    # Not a boolean.
                    error_msg = "Expected a boolean value."
                    self.inline_errors[path] = error_msg

        # Validate strings and other values.
        def _validate_other(self, values):
            for value in values:
                if type(value).__name__ != data_type:
                    # Other type mismatch.
                    error_msg = f"Expected type '{data_type}' instead of '{type(value).__name__}'"
                    self.inline_errors[path] = error_msg

        #
        #

        # Skip validation if no schema is provided.
        if not self.schema:
            return

        # We're validating he currently selected
        # value, unless a path is provided.
        path = path if path else self.path

        # Check if the key is in the schema.
        # This is only relevant when you are loading a file.
        if not self.value_from_path(data=self.schema, path=path):
            self.illegal_errors[path] = 1
            return

        # Get requirements from schema.
        required = self.value_from_path(
            data=self.schema, path=path + ('required',))
        data_type = self.value_from_path(
            data=self.schema, path=path + ('type',))

        # Get the value from the data.
        value = self.value_from_path(path=path)

        # If the data type is a list, we need
        # to validate by the type in the list.
        is_list = False
        if type(data_type) == list:
            data_type = data_type[0]
            values = value
            is_list = True
        else:
            values = [value]

        # Validate lists.
        if is_list and values != '' and type(values) != list:
            error_msg = "Expected a list."
            self.inline_errors[path] = error_msg
            return

        if len(values) == 1 and values[0] == '' and required:
            # Required field is empty.
            error_msg = "Value required."
            self.inline_errors[path] = error_msg
        elif data_type == 'enum':
            # Validate enums.
            _validate_enum(self, values)
        elif data_type == 'int' or data_type == 'float':
            # Validate numbers.
            _validate_number(self, values, data_type)
        elif data_type == 'bool':
            # Validate booleans.
            _validate_bool(self, values)
        else:
            # Validate all other values (str).
            _validate_other(self, values)

    # Remove all illegal fields from self.data.
    def format_file(self):
        for path in self.illegal_errors:
            self.remove_key_from_path(path=path)
        self.illegal_errors = {}

    # Write changes to JSON file.
    def update_file(self):
        if bool(self.inline_errors):
            self.error = 'Some of your changes are invalid.'
        else:
            try:
                # raise ValueError("Invalid input.") # To simulate error.
                if self.demo:
                    return self.term.green('Changes saved to JSON file.') + self.term.yellow('\nNote: your were in demo mode so your changes have not really been saved.')
                else:
                    print(self.json_file_path)
                    with open(self.json_file_path, 'w', encoding='UTF-8') as f:
                        print(1)
                        json.dump(self.data, f, indent=4)
                        print(2)
                        f.close()
                        print(3)
                    return self.term.green('Changes saved to JSON file.')
            except ValueError as error_msg:
                self.error = 'There was an error saving your changes.'
                print(error_msg)


if __name__ == '__main__':
    ej = EditJson()
    ej('.demo/demo.json', '.demo/demo-schema.json')
