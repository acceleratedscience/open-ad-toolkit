import ast
import os
import re
import operator
import json
import math
import copy
from functools import reduce
from threading import Timer
from collections import deque, OrderedDict

# Core
import jsonschema
from blessed import Terminal

# Importing our own plugins.
# - - -
# This is a temporary solution until every plugin is
# available as a public pypi package.
# - - -
# In production, every plugin has its own package directory
# with __init__.py included, which all live inside the plugin
# directory. The plugin directory has an __init__.py itself
# which lets our plugins import eachother via the parent package
# using ..cousin_module.
# In development mode, each plugin directore lives in its
# develoment directory, which are meant to be sibling dirs so
# they can import one another using sys.path.append('../')
try:
    # Used in production setup.
    from ..style_parser import print_s, style, strip_tags, a_len, a_textwrap
except BaseException:
    # Used in isolated plugin development setup.
    import sys

    sys.path.append("../")
    from style_parser import print_s, style, strip_tags, a_len, a_textwrap

# Debug modes:
# 1: Display logger
# 10: Display logger in overlay
# 2: Display line numbers (zero based)
# 3: Display content/terminal height
# 4: Display data
# 5: Log key codes (combine with 1)
# 6: Display navigational parameters
DEBUG = []


class EditJson:
    """
    Edit a JSON file within the CLI.
    --------------------------------
    Author: Moenen Erbuer - moenen.erbuer@ibm.com
    v0.0.0-beta5 / Last update: Sep 29, 2023

    Description:
        This module lets you edit JSON from within the CLI with an
        intuitive GUI-like interface. When a json schema is included,
        we'll use it to validate input and to display field descriptions.

    Some useful links about json schema:
        - Overview: https://jsonschema.org
        - Data type reference: https://json-schema.org/understanding-json-schema/reference/index.html
        - Autogenerate your schemas: https://www.jsonschema.net/app/schemas/0


    Parameters:
        json:str, default=None
            Path to a JSON file. If no path is provided, we will ask for one.
        schema:str, default=None
            Path to a json schema (see jsonschema.org) that contains data types
            and description per field. Used for validation & instructions.
        template:str, default=None
            When set and no JSON file was found, a new file will be created with
            the template's content.
        title:str, default=None
            Custom title on top of the JSON editor.
        new:bool, default=False
            When set to True and no JSON file was found, a new file will be created
            from the schema.
        strict:bool, default=False
            When set to True, a warning message will be displayed on top whenever
            any schema fields are missing, regardless of whether they're required.
            By default, the warning message will only show whenever there's any
            required fields missing.
        read_only:bool, default=False
            Opens the file in read-only mode.
        demo:bool, default=False
            Used in isolation of any other parameter, this launches the JSON editor
            in demo mode, with sample data showcasing capabilities.



    Usage:
        `
        import edit_json

        # Let user select a JSON file.
        edit_json()

        # Load JSON file without schema.
        edit_json('_sample_files/sample.json')

        # Load JSON file with a schema. Throw error when no JSON file is found.
        edit_json('_sample_files/sample.json', '_sample_files/sample-schema.json')

        # Create new JSON file from the schema when none is found.
        edit_json('_sample_files/sample.json', '_sample_files/sample-schema.json', new=True)

        # Create new JSON file from the template when none is found.
        edit_json('_sample_files/sample.json', template='_sample_files/sample-template.json')

        # Create new JSON file from the schema after asking user to select a path
        # for it to be saved. Not that file creation only happens on exit.
        edit_json(schema='_sample_files/sample-schema.json', new=True)

        # Strict mode: Missing non-required fields will also trigger
        # a warning, including instructions on how to format the file.
        # edit_json('_sample_files/sample.json', '_sample_files/sample-schema.json', strict=True)

        # Run the editor in demo mode with dummy data.
        edit_json(demo=True)

        # Consult the docstring for usage info.
        help(edit_json)
        `
    """

    screens = {}
    active_screen = None  # This lets us render different screens, like help, save, etc.

    def __init__(self):
        # Screen management (help, confirm_exit, maybe more later)
        self.active_screen = None  # To store an active screen's class.
        screen_classes = [ConfirmExitScreen, AllCommandsScreen]
        for ScreenClass in screen_classes:
            self.screens[ScreenClass.name] = ScreenClass(self)

        # Terminal
        self.term = Terminal()
        try:
            term_size = os.get_terminal_size()
            self.width = term_size.columns
            self.height = term_size.lines - 1  # TODO: Figure out why we need -1 here, blessed not using full height?
        except:
            self.width = 80
            self.height = 24
            pass  # is not supporting get terminal

        self.display_height = 0  # The height of the output

        # Files
        self.json_path = ""  # Path to JSON file
        self.json_filename = ""  # Name of JSON file
        self.schema_path = ""  # Path to schema file

        # Data storage
        self.text_input = ""  # To store typing until you press enter.
        self.data = {}  # To store JSON data until it gets saved to file.
        self.addresses = []  # List of addresss tuples - see store_addresses()
        self.selected_address = None  # The address (tuple with key hierarchy) for the selected field.
        self.schema = None  # Dictionary with schema information.

        # Editor UI
        self.edit_mode = False
        self.cursor_pos = 0
        self.instructions_main = ""
        self.empty_placeholder = "<soft>-</soft>"
        self.blank = False  # Indicates if JSON file is generated (and thus file not yet saved)
        self.blank_from_schema = False  # Indicates if blank JSON is generated from schema.
        self.blank_from_template = False  # Indicates if blank JSON is generated from template.
        self.ui_width = 0  # The width of the instructions and separator lines on top, used to wrap text.
        self.left_col_width = 0  # Width of the key column.
        self.right_col_width = 0  # Width of the value column.
        self.indent_str = "<bright_black>|   </bright_black>"  # Indentation per level.
        self.indent_str_len = 4
        self.skip_lines = 0  # Used to emulate viewport scrolling.
        self.read_only = False  # Block editing
        self.schema_mode = False  # When displaying the schema
        self.show_path = False  # Display the JSON file path
        self.sep = style("<soft>" + "-" * self.width + "</soft>", nowrap=True)

        # Errors
        self.error = ""  # General error message on top.
        self.inline_errors = {}  # Inline error messages (invalic input according to schema).
        self.illegal_errors = {}  # Illegal field flags (field not present in schema).

        # Navigation
        self.header_height = 4
        self.first_line = 0  # Position of first line of data, below header.
        self.line_nr = 0  # Used to render the UI line by line.
        self.selected_index = 0  # Index of field in focus.
        self.selected_line_nr = None  # The line number of the field in focus.
        self.index_to_line_nr_map = {}  # Used to deduct the line number from the selected index.
        self.selected_key = ""  # The JSON key in focus.
        self.selected_value = ""  # The JSON value in focus.
        self.selected_help = ""  # The help string for the key in focus.
        self.show_help = False  # Toggle the help string for the key in focus.
        self.listen_for_action = False  # When you type `:` we wait for `s/x/f...` commands.
        self.esc_code = False  # Used to detect alt-left/right key.

        # Options
        self.template = None  # A template JSON file we'll copy all data from.
        self.title = None  # Displays on top of header, replacing default "You are editing xyz.json".
        self.new = False  # Create new JSON file from schema when if it doesn't exist.
        self.strict = None  # When True, an error is displayed when the JSON doesn't fully match the schema.
        self.read_only = False  # Open file in read-only mode.
        self.demo = False  # Demo mode - uses sample JSON file and blocks file saving.
        #
        self.options_hist = []  # When jumping from file to file, this lets you go back.

        # Dev
        self.log_values = []  # Used by self.log()
        self.og_display_lines = []  # Used to store display lines before self.log() manipulates them.

    def __call__(self, *args, **kwargs):
        success = self.open_json_file(*args, **kwargs)
        if success:
            self.render_ui()

    ##
    # Initiation
    ##

    # Open JSON file and store all options.
    def open_json_file(
        self,
        json=None,  # Required: JSON file path.
        schema=None,  # Optional: schema for the JSON file.
        # Options - see above.
        template=None,
        title=None,
        new=False,
        strict=False,
        read_only=False,
        demo=False,
    ):
        self.json_path = json
        self.schema_path = schema
        template_path = template
        self.title = title
        self.new = new
        self.strict = strict
        self.read_only = read_only
        #
        # Store all options in a list, so when you jump to another
        # JSON file with different options, you can easily go back.
        self.options_hist.append([json, schema, template, title, new, read_only, demo])  # Keep in sync with arguments!

        # Reset
        self.selected_index = 0

        # Load demo mode data.
        if demo:
            dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "demo")
            self.json_path = os.path.join(dir, "demo.json")
            self.schema_path = os.path.join(dir, "demo-schema.json")
            self.demo = True
            self.title = "You are in demo mode."

        # Load schema file.
        if self.schema_path:
            self.schema = self.load_schema_file(self.schema_path)
            if not self.schema:
                return
        elif self.new is True:
            print_s("<error>A schema is required when new=True.</error>", nowrap=True)
            return
        else:
            self.schema = None

        # Load template file.
        if template_path:
            template_success = self.load_template_file(template_path)
            if not template_success:
                return

        # Request JSON path is none provided
        self.json_path = self.json_path if self.json_path else self.request_path()

        # Validate & load JSON data.
        data, file_exists, file_valid = self.load_json_file(self.json_path)

        # Create blank JSON file if need be.
        if not file_exists:
            if self.new:
                data = self.create_file_from_schema()
                self.blank = True
                self.blank_from_schema = True
            elif template_path:
                data = self.create_file_from_template()
                self.blank = True
                self.blank_from_template = True
            else:
                return
        elif not file_valid:
            return

        # Store data.
        self.data = data
        self.json_path = self.json_path
        self.json_filename = self.filename_from_path(self.json_path)

        # Parse JSON data to set up UI.
        self.store_addresses()
        self.set_left_col_width()
        self.validate()

        # For debugging.
        # Print data & key addresses (replaces UI).
        # - - - - - - - - - - - - - -
        if 4 in DEBUG:
            print("\ndata (top level):\n-----\n" + "\n".join(self.data))
            print("\nAddresses:\n-------")
            for self.json_path in self.addresses:
                print(self.json_path)
            print("\n")
            return
        # - - - - - - - - - - - - - -
        return True

    # Ask user for path to JSON file.
    def request_path(self):
        while True:
            json_path = input(style(f"<yellow>Enter JSON file path: </yellow>")).strip()
            if json_path:
                break
        return json_path

    # Load schema file.
    def load_schema_file(self, schema_path):
        schema = None
        try:
            with open(schema_path, encoding="UTF-8") as file:
                schema = json.load(file)
                # raise FileNotFoundError("Bla bla bla.")  # For testing
        except FileNotFoundError as err:
            print_s(
                f"<error>Schema file not found: <yellow>{schema_path}</yellow></error>\n<soft>{err}</soft>", nowrap=True
            )
        except json.JSONDecodeError as err:
            print_s(
                f"<error>Invalid schema file: <yellow>{schema_path}</yellow></error>\n<soft>{err}</soft>", nowrap=True
            )

        if schema:
            # Expand schema $refs.
            schema = self.expand_schema_refs(schema)

            # Validate schema.
            try:
                jsonschema.validate(schema=schema, instance={})
            except jsonschema.exceptions.SchemaError as err:
                print_s(
                    f"<error>Invalid schema</error> <yellow>{schema_path}</yellow>\n<soft>{err.message}</soft>",
                    nowrap=True,
                )
                schema = None
            except jsonschema.exceptions.ValidationError as err:
                pass

        return schema

    # A key's schema can be a reference to another schema.
    # Whenever we're looking at a key's schema this function
    # needs to be invoked to make sure that remote schemas
    # are looked up and expanded.
    # https://json-schema.org/understanding-json-schema/structuring.html#ref
    def expand_schema_refs(self, schema):
        def _recursion(schema_props):
            for key, prop in schema_props.items():
                if "$ref" in prop:
                    # self.log('\n//', key)
                    # Pull up subschema from $defs.
                    if prop["$ref"][:8] == "#/$defs/":
                        ref_key = prop["$ref"].split("/")[-1]
                        if "$defs" in schema:
                            schema_props[key] = schema["$defs"][ref_key]
                            # self.log('A >>>', schema['$defs'][ref_key])

                    # Pull up schema from path.
                    elif prop["$ref"][0] == "/":
                        path = prop["$ref"]
                        result = {}

                        def _recursion_find_schema_ref(schema_props2):
                            nonlocal result
                            for i, (key, prop2) in enumerate(schema_props2.items()):
                                # self.log('@', i, key)
                                if "properties" in prop2 and isinstance(prop2["properties"], dict):
                                    # self.log('-->')
                                    abort = _recursion_find_schema_ref(prop2["properties"])
                                    if abort:
                                        break
                                elif "$id" in prop2 and prop2["$id"] == path:
                                    prop2_clone = copy.copy(prop2)
                                    del prop2_clone["$id"]
                                    result = prop2_clone
                                    # self.log('B >>>', result)
                                    return result
                            return result

                        schema_props[key] = _recursion_find_schema_ref(schema["properties"])

                    # Pull up recursive schema. TODO!
                    # https://json-schema.org/understanding-json-schema/structuring.html#extending-recursive-schemas
                    elif prop["$ref"] == "#":
                        # schema_prop = schema
                        pass

        _recursion(schema["properties"])
        return schema

    # Load JSON file.
    def load_json_file(self, json_path):
        data = None
        file_exists = False
        file_valid = True
        try:
            with open(json_path, "r", encoding="UTF-8") as file:
                data = json.load(file)
                file_exists = True
        except FileNotFoundError:
            if not self.new and not self.template:
                print_s(f"<error>JSON file not found: <yellow>{json_path}</yellow></error>", nowrap=True)
        except json.JSONDecodeError as err:
            file_exists = True
            file_valid = False
            print_s(f"<error>Invalid JSON file: <yellow>{json_path}</yellow></error>\n<soft>{err}</soft>", nowrap=True)
        except BaseException as err:
            print_s(f"<error>Somethig went wrong</error>\n<soft>{err}</soft>", nowrap=True)
        return data, file_exists, file_valid

    # Load template file.
    def load_template_file(self, template_path):
        success = False
        try:
            with open(template_path, encoding="UTF-8") as file:
                # raise FileNotFoundError("Bla bla bla.")  # For testing
                self.template = json.load(file)
                success = True
        except FileNotFoundError as err:
            print_s(
                f"<error>Template file not found: <yellow>{template_path}</yellow></error>\n<soft>{err}</soft>",
                nowrap=True,
            )
        except json.JSONDecodeError as err:
            print_s(
                f"<error>Invalid template file: <yellow>{template_path}</yellow></error>\n<soft>{err}</soft>",
                nowrap=True,
            )
        return success

    # Create JSON file from template.
    def create_file_from_template(self):
        data = self.template
        return data

    # Create JSON file from schema.
    def create_file_from_schema(self):
        def _recursion_create_one_level(schema, data):
            if "properties" in schema:
                for key, val in schema["properties"].items():
                    if "properties" in val:
                        data[key] = {}
                        _recursion_create_one_level(val, data[key])
                    else:
                        data[key] = val["default"] if "default" in val else ""

        #
        #

        data = {}
        _recursion_create_one_level(self.schema, data)
        return data

    # Split file name from file path
    def filename_from_path(self, path):
        path_tree = path.split("/")
        return path_tree[len(path_tree) - 1]

    # Store key addresses.
    # - - -
    # Generated flat array with a tuple per key,
    # describing its address in the disctionary.
    # This is used to navigate between fields when
    # pressing UP/DOWN keys.
    # - - -
    # For example:
    # data = { foo: 1, bar: { a: 2, b: 3 } }
    # --> [('foo'), ('bar', 'a'), ('bar', 'b')]
    def store_addresses(self):
        self.addresses = []

        def _recursion(data, parent_keys=tuple()):
            for key, val in data.items():
                is_dict = isinstance(val, dict)
                if is_dict:
                    _recursion(val, parent_keys=parent_keys + (key,))
                else:
                    self.addresses.append(parent_keys + (key,))

        _recursion(self.data)

    # Measure the length of all keys on all levels,
    # accounting for indentation, and set the left
    # column width so it fits everything, but cap
    # it at 50 characters, or half the terminal width
    # if that's less.
    def set_left_col_width(self):
        def _recursion(data, level=0):
            max_key_length = 0
            for i, (key, value) in enumerate(data.items()):
                is_dict = isinstance(value, dict)
                if is_dict:
                    max_key_length = max(
                        _recursion(value, level=level + 1), len(key) + self.indent_str_len * level, max_key_length
                    )
                else:
                    max_key_length = max(len(key) + self.indent_str_len * level, max_key_length)
            return max_key_length

        max_key_length = _recursion(self.data)
        upper_limit = min(math.floor(self.width / 2), 50)
        self.left_col_width = min(max_key_length + 2, upper_limit)
        self.right_col_width = self.width - self.left_col_width - 2  # Two chars used for text field padding.

    ##
    # Validation
    ##

    # Validate all fields according to the schema.
    def validate(self):
        if not self.schema:
            return
        validator = jsonschema.Draft202012Validator(self.schema)
        errors = validator.iter_errors(self.data)
        for i, err in enumerate(errors):
            # When an array's value doesn't match the schema,
            # the err.path will return the index of the faulty
            # value as last item in the path list. Because we
            # only display one error per field, this value
            # needs to be stripped from the path.
            error_address = list(err.path)
            is_array_value_error = len(err.path) and isinstance(err.path[-1], int)
            if is_array_value_error:
                error_address.pop()
            error_address = tuple(error_address)
            self.inline_errors[error_address] = err.message

        # Flag required fields with empty values.
        required_errors = self.validate_required(self.data, self.schema)
        self.log(required_errors)
        for err in required_errors:
            error_address = tuple(err["path"])
            self.inline_errors[error_address] = err["message"]

        # Detect and store any illegal keys (i.e. keys not present in schema).
        for i, address in enumerate(self.addresses):
            schema_address = self.address_to_schema_address(address)
            is_illegal = not bool(self.value_from_address(self.schema, schema_address))
            if is_illegal:
                self.illegal_errors[address] = 1
                del self.addresses[i]

    # According to the jsconschema spec, an empty required field validates
    # as True, because the key is present. We want to validate it as False.
    # https://stackoverflow.com/a/54426110
    def validate_required(self, data, schema, address=()):  # TODO: Make this recursive
        errors = []

        if isinstance(data, dict):
            for required_field in schema.get("required", []):
                if required_field in data and data[required_field] == "":
                    req_address = address + (required_field,)
                    errors.append(
                        {
                            "message": f"'{required_field}' is a required property",
                            "path": req_address,
                        }
                    )

            for key, value in data.items():
                key_address = address + (key,)
                # self.log(key_address)
                if isinstance(value, dict) and "properties" in schema:
                    if key in schema["properties"]:
                        errors.extend(self.validate_required(value, schema["properties"][key], key_address))
                else:
                    errors.extend(self.validate_required(value, schema, key_address))

        # # This is for when your data is a list of objects instead of an object... not sure if we need,
        # and if we do this needs to be tested.
        # elif isinstance(data, list):
        #     for idx, item in enumerate(data):
        #         # new_path = f"{path}[{idx}]"
        #         new_path = path + (idx,)
        #         print('#', new_path)
        #         errors.extend(validate_required(item, schema, new_path))

        return errors

    ##
    # Rendering
    ##

    # Render the UI.
    def render_ui(self):
        # Display instructions on top.
        def get_instructions_main():
            output = [
                "<reverse> ▲ </reverse> <reverse> ▼ </reverse> Navigate",
                "<reverse> ▶ </reverse> Property info" if self.schema else None,
                "<reverse> ESC </reverse> Exit",
                "<cmd>:ENTER</cmd> Save",
                "<cmd>:?</cmd> Help",
            ]
            output = [item for item in output if item]
            output = "   ".join(output)
            if not self.ui_width:
                # Set ui width to width of the instruction string, or at least 80 characters.
                self.ui_width = max(min(len(strip_tags(output)), self.width), 80)
            return style(output, nowrap=True)

        def get_instructions_read_only():
            output = [
                "<reverse> ESC </reverse> " + ("Go back" if self.schema_mode else "Exit"),
            ]
            output = "   ".join(output)
            return style(output, nowrap=True)

        self.instructions_main = get_instructions_main()
        self.instructions_read_only = get_instructions_read_only()
        exit_msg = None

        with self.term.fullscreen(), self.term.cbreak(), self.term.hidden_cursor():
            while True:
                self.render_one_frame()

                # Process user input (ks = keystroke)
                # Note: the while loop waits for the inkey input to be released.
                try:
                    # raise BaseException('Lorem Ipsum')
                    ks = self.term.inkey(timeout=3, esc_delay=0)
                    exit_msg = self.process_key_stroke(ks)
                    if exit_msg:
                        break
                except KeyboardInterrupt:
                    exit_msg = self.confirm_exit("discard")
                    if exit_msg:
                        break
                # except BaseException as err:
                #     exit_msg = style(f'<error>Something went wrong</error>\n<soft>{err}</soft>')
                #     break

        # Print message after you exited: "Changes saved/discarded/none"
        if exit_msg and exit_msg != True:
            print_s(exit_msg, nowrap=True)

    # Render single frame every time a key is pressed (see inkey).
    def render_one_frame(self):
        display = []
        self.first_line = 0

        # We count the lines of everything added to display
        # so we can keep track of the line numbers of our keys.
        # Note: render_list() has its own version of this function.
        # This one is used for general non-list content like titles
        # and general error messages.
        def display_append(output, full_width=False):
            if not output:
                return

            width = self.width if full_width else self.ui_width
            lines = a_textwrap(output, width=width, drop_whitespace=False).splitlines()
            for line in lines:
                display.append(line)
                self.first_line += 1

        # Store selected value.
        if not self.read_only:
            try:
                # Re-enable the lines below
                self.selected_address = self.addresses[self.selected_index]
                self.selected_key = self.selected_address[-1]
                self.selected_value = self.value_from_address(self.data, self.selected_address)
            except BaseException as err:
                print_s(f"<error>Something went wrong rendering the UI</error>\n<soft>{err}</soft>", nowrap=True)
                return

        # Render help string for selected field.
        if self.schema and not self.read_only:
            selected_schema_address = self.address_to_schema_address(self.selected_address)
            schema_prop = self.value_from_address(self.schema, selected_schema_address)

            if not schema_prop or not isinstance(schema_prop, dict):
                self.selected_help = None
            else:
                help_str = schema_prop["description"] if schema_prop and "description" in schema_prop else ""
                type_str = schema_prop["type"] if schema_prop and "type" in schema_prop else ""
                if isinstance(type_str, list):
                    type_str = "/".join(type_str)
                range_str = ""
                multiple_of_str = ""

                # Regarding mathematical annotating of inclusive & exclusive
                # ranges with parentheses or brackets - [0,10) / (0,10] / etc:
                # https://stackoverflow.com/a/4396303
                if type_str == "integer" or type_str == "number":
                    # Range
                    if "minimum" in schema_prop:
                        if "maximum" in schema_prop:
                            range_str = f" [{schema_prop['minimum']}-{schema_prop['maximum']}]"
                        elif "exclusiveMaximum" in schema_prop:
                            range_str = f" [{schema_prop['minimum']}-{schema_prop['exclusiveMaximum']})"
                        else:
                            range_str = f">={schema_prop['minimum']}"
                    elif "exclusiveMinimum" in schema_prop:
                        if "maximum" in schema_prop:
                            range_str = f" ({schema_prop['exclusiveMinimum']}-{schema_prop['maximum']}]"
                        elif "exclusiveMaximum" in schema_prop:
                            range_str = f" ({schema_prop['exclusiveMinimum']}-{schema_prop['exclusiveMaximum']})"
                        else:
                            range_str = f">{schema_prop['exclusiveMinimum']}"
                    elif "exclusiveMaximum" in schema_prop:
                        range_str = f"<{schema_prop['exclusiveMaximum']}"
                    elif "maximum" in schema_prop:
                        range_str = f"<={schema_prop['maximum']}"

                    # Multiple of
                    if "multipleOf" in schema_prop:
                        multiple_of_str = f":{schema_prop['multipleOf']}"
                self.selected_help = (
                    f"<green>{type_str.upper()}{multiple_of_str}{range_str}</green><yellow/> {help_str}"
                    if type_str and help_str
                    else f"<green>{type_str.upper()}</green>"
                    if type_str
                    else str(help_str)
                )  # noqa

        if not self.active_screen or self.active_screen.keep_title:
            # Generate status flag.
            if self.read_only:
                status_flag = "<on_red> READ ONLY </on_red> "
            else:
                status_flag = "<on_green> NEW FILE </on_green> " if self.blank else f"<on_yellow> EDITING </on_yellow> "

            # Add status flag + Title
            if self.schema_mode:
                title_output = style(f"{status_flag}<yellow>You're viewing the schema</yellow>", nowrap=True)
            elif self.title:
                title_output = style(f"{status_flag}<yellow>{self.title}</yellow>", nowrap=True)
            else:
                title_output = style(f"{status_flag}<yellow>{self.json_filename}</yellow>", nowrap=True)
            display_append(title_output)

            # ------------
            display_append(self.sep, True)

        # DISPLAY ACTIVE SCREEN (confirm exit, help...)
        if self.active_screen:
            display.append(self.render_active_screen())

        # DISPLAY DATA
        else:
            # Instructions
            if self.listen_for_action:
                instructions_str = ":"
            elif self.read_only:
                instructions_str = self.instructions_read_only
            else:
                instructions_str = self.instructions_main
            display_append(instructions_str, True)

            # ------------
            display_append(self.sep, True)

            # JSON file path
            if self.show_path:
                path_output = style(f"\U0001F4C1 File path: <magenta>{self.json_path}</magenta>", nowrap=True)
                if 1:
                    path_output += style(" / Type <cmd>:o</cmd> to open", nowrap=True)
                display_append(path_output)
                display_append(self.sep, True)

            # General error - invalid value (after saving attempt)
            error_general = None

            # Detect discrepancies with schema.
            missing_keys, illegal_keys = self.detect_discrepancies()
            if missing_keys or illegal_keys:
                error_general = style(
                    f"<error>This JSON is not matching the schema. Type <cmd>:f</cmd> to reformat it.</error>",
                    nowrap=True,
                )
            elif self.error:
                error_general = style(f"<error>{self.error}</error>")
            # General error - illegal values
            elif self.illegal_errors:
                error_general = style(
                    "<error>This JSON file does not match the provided schema.\nType <cmd>:c</cmd> to clear out illegal values or <cmd>:f</cmd> to fully reformat the file.</error>",
                    nowrap=True,
                )  # noqa

            if error_general:
                display_append(error_general)
                display_append(self.sep, True)

            # The line we're currently rendering.
            self.line_nr = self.first_line

            # The line that's selected.
            if self.selected_line_nr is None:
                self.selected_line_nr = self.line_nr

            # Print list
            display += self.render_list()

        # ------------
        display_append(self.sep, True)

        # For debugging
        # Log data below list using self.log()
        # - - - - - - - - - - - - - -
        if 1 in DEBUG:
            display.append(f"\n----")
            for val in self.log_values:
                display.append(val)  # Debug
        # - - - - - - - - - - - - - -

        # Auto-scroll behavior.
        # Trim the display output to fit within the console viewport.
        initial_display_len = len(display)
        if self.read_only:
            first_line = self.skip_lines
            last_line = self.skip_lines + self.height
            display = display[first_line:last_line]
        else:
            if len(display) > self.height:
                first_line = self.skip_lines
                last_line = self.skip_lines + self.height
                if self.selected_line_nr >= last_line - 4:
                    max_skip = initial_display_len - self.height
                    self.skip_lines = min(self.selected_line_nr - self.height + 1 + 4, max_skip)
                elif self.selected_line_nr <= first_line + self.first_line:
                    self.skip_lines = max(self.selected_line_nr - self.first_line, 0)
                first_line = self.skip_lines
                last_line = self.skip_lines + self.height
                display = display[first_line:last_line]

        # Clear the console & render display.
        print(self.term.home + self.term.clear)

        # Print self.log as overlay
        # INSTABLE IMPLEMENTATION - TO BE CLEANED UP
        if 10 in DEBUG:
            for i, val in enumerate(self.log_values):
                val = style(f"<reverse> {val} </reverse>", nowrap=True)
                # val = style(f'{"":>{self.width - a_len(val) - 2}}<reverse> {val} </reverse>', nowrap=True) ???
                if len(display) > i:
                    display[i] = val

        print("\n".join(display))
        self.display_height = len(display)

        # For debugging.
        # - - - - - - - - - - - - - -
        if 3 in DEBUG:
            # Display output height vs terminal height.
            print(f"\nOutput/Terminal: {initial_display_len}/{self.height}")
        if 6 in DEBUG:
            # Display navigational data.
            print("\nEdit mode:", self.edit_mode)
            print("Cursor:", self.cursor_pos)
            print("Listen:", self.listen_for_action)
            print("Input:", self.text_input)
            print("Index:", self.selected_index)
            print("Path:", self.selected_address)
            print("Selected key:", self.selected_key)
            print("Selected value:", self.selected_value)
            print("Selected help:", self.selected_help)
            print("Inline errors: ", self.inline_errors)
            print("\n")
        # - - - - - - - - - - - - - -

    # Render interactive list of key-value pairs.
    def render_list(self):
        display = []
        index = 0

        # Sometimes the > arrow wraps to new line,
        # in which case the help string shouldn't
        # start with a line break.
        help_str_new_line = True

        def _render_one_line(data, indent_level=0, address=()):
            nonlocal index
            for i, (key, value) in enumerate(data.items()):
                value = self.empty_placeholder if not value and value != 0 else value
                is_dict = isinstance(value, dict)
                indent = self.indent_str * indent_level
                indent_len = self.indent_str_len * indent_level
                key_address = address + (key,)
                self.index_to_line_nr_map[index] = self.line_nr
                is_illegal = key_address in self.illegal_errors

                # For debugging
                # Display index / line numbers / line nr in focus.
                # - - - - - - - - - - - - - -
                line_nr_debug = ""
                line_nr_debug_blank = ""
                right_col_width = self.right_col_width
                if 2 in DEBUG:
                    index_str = "--" if is_dict or is_illegal else f"{index:02}"
                    line_nr_debug = f"{index_str}:{self.line_nr:02}:{self.selected_line_nr:02} "
                    line_nr_debug_blank = " " * len(line_nr_debug)
                    right_col_width = self.right_col_width - len(line_nr_debug)
                # - - - - - - - - - - - - - -

                # Parent field (no value).
                if is_dict:
                    display.append(
                        style(line_nr_debug + indent + f"<soft><underline>{key}:</underline></soft>", nowrap=True)
                    )
                    self.line_nr = self.line_nr + 1
                    _render_one_line(value, indent_level=indent_level + 1, address=key_address)

                # Field with value.
                else:
                    line_gap = (
                        self.left_col_width - self.indent_str_len * indent_level
                    )  # The length of the gap between the indentation and the value. noqa
                    is_selected = self.selected_address == key_address
                    illegal_flag = style("<error_reverse> ILLEGAL </error_reverse>") if is_illegal else ""
                    left_col_filler = (
                        line_nr_debug_blank + style(indent, nowrap=True) + " " * (self.left_col_width - indent_len)
                    )

                    # Selected row.
                    if is_selected and not illegal_flag:
                        key_str = line_nr_debug + indent + f"{(key + ':'):<{line_gap}}"

                        # Edit mode.
                        if self.edit_mode:
                            value_str = __render_edit_mode(left_col_filler, right_col_width)

                        # Editable.
                        else:
                            value_str = __render_editable(value, left_col_filler, right_col_width)

                        # Helper text.
                        if self.show_help and self.selected_help:
                            value_str += __render_helper_text(left_col_filler, right_col_width)

                        # Store line ouput.
                        line_output = style(key_str + value_str, nowrap=True)

                        # In case the key is wider than the left column, we trim it to fit in edit mode.
                        key_str_len = a_len(style(key_str))
                        if key_str_len > self.left_col_width:
                            trimmed_key = key[: self.left_col_width - indent_len - 3] + "…"
                            key_str = line_nr_debug + indent + f"{(trimmed_key + ':'):<{line_gap}}"
                            line_output = style(key_str + value_str, nowrap=True)

                        __display_append(line_output)

                    # Non-selected row.
                    else:
                        # Key
                        key_str = line_nr_debug + indent + f"<soft>{(key + ':'):<{line_gap}}</soft> "

                        # Value
                        if illegal_flag:
                            # Illegal fields - [ILLEGAL] + grey strikethrough text
                            value = "" if value == self.empty_placeholder else value
                            value_str = __render_illegal(illegal_flag, value, left_col_filler, right_col_width)
                            key_value_space_len = len(strip_tags(key_str)) - len(strip_tags(indent)) - len(key) - 1
                            if key_value_space_len > 1:
                                key_str = key_str[:-1]
                        elif self.inline_errors and self.addresses[index] in self.inline_errors:
                            # Invalid value - red
                            value_str = __render_invalid(value, left_col_filler, right_col_width)
                        else:
                            # Regular field.
                            value_str = __render_regular(value, left_col_filler, right_col_width)

                        # Store line ouput.
                        line_output = style(key_str + value_str, nowrap=True)

                        # In case the key + value don't fit on one line, we trim the key.
                        line_output_len = a_len(line_output.splitlines()[0])
                        if line_output_len > self.width:
                            trim = self.width - line_output_len - 2
                            trimmed_key = key[:trim] + "…"
                            key_str = line_nr_debug + indent + f"<soft>{(trimmed_key + ':'):<{line_gap}}</soft> "
                            line_output = style(key_str + value_str, nowrap=True)

                        __display_append(line_output)

                    # Inline error message.
                    if not illegal_flag:
                        if self.inline_errors and self.addresses[index] in self.inline_errors:
                            line_output = __render_error_msg(left_col_filler, right_col_width)
                            __display_append(line_output)

                        index += 1

        #
        #

        def __render_edit_mode(left_col_filler, right_col_width):
            value_str = a_textwrap(self.text_input, width=right_col_width, drop_whitespace=False)  # noqa
            lines = value_str.splitlines()
            chars_before_cursor = 0

            # AVERT YOUR EYES - THIS IS A HACK
            # - - - - - - - - - - - - - - - - - - - - - - - - - -
            # In a regular text editor, you will never have a line
            # that starts with a space, even when the previous line
            # ends with one or more spaces that should be wrapped.
            # This is done by using zero-width characters, however
            # this doesn't work in the CLI because we can't reliably
            # print the standard zero-width character (\u200b).
            # I was able to achieve it by using the antiquated and
            # thus non-reliable PADDING character (\x80), but then
            # realized that using a zero-width character also
            # complicates things by making it harder to calculate
            # filler_right correctly. We can't simply trim the
            # whitespace either because that messes with the cursor
            # position. The easy way out would be to simply allow
            # new lines to start with a space but this is a really
            # poor look. So instead we're checking every line, and
            # in case there is ONE space that can be trimmed, we
            # add it to the previous line, and remove the one-char
            # right textbox padding for that line. Unfortunately
            # this doesn't work when a line is starting with more
            # than one space, so if that's the case we don't have a
            # choice but to just leave it. This is acceptable though
            # because it is purely aesthetical and our solution
            # covers 99% of cases.
            pad_right = []
            for i, line in enumerate(lines):
                pad_right += " "
                line_stripped = line.lstrip()
                stripped_len = len(line) - len(line_stripped)
                if stripped_len == 1:
                    lines[i] = line_stripped
                    pad_right[i - 1] = ""
                    lines[i - 1] += " "

            for i, line in enumerate(lines):
                is_last_line = i == len(lines) - 1
                if self.cursor_pos - chars_before_cursor < 0:
                    # Lines after the line with cursor
                    pass
                elif self.cursor_pos - chars_before_cursor >= len(line) and not is_last_line:
                    # Lines before the line with cursor.
                    chars_before_cursor += len(lines[i])
                else:
                    # Line with cursor.
                    pos = self.cursor_pos - chars_before_cursor
                    filler_right = right_col_width - len(line)
                    cursor_str = line[pos : pos + 1]

                    # When the cursor is at the very end we fill
                    # it up with a space so you can still see it.
                    if len(cursor_str) == 0:
                        cursor_str = " "
                        filler_right -= 1

                    lines[i] = (
                        line[:pos] + f"<on_yellow>{cursor_str}</on_yellow>" + line[pos + 1 :] + filler_right * " "
                    )
                    chars_before_cursor += len(line)
            value_str = (
                "<edit> "
                + f"</edit>\n{left_col_filler}<edit> ".join(
                    [f"{line:<{right_col_width}}{pad_right[i]}" for i, line in enumerate(lines)]
                )
                + "</edit>"
            )  # noqa
            return value_str

        def __render_editable(value, left_col_filler, right_col_width):
            nonlocal help_str_new_line
            show_help_str = (
                ""
                if not self.selected_help
                else style(" <yellow>▼</yellow>")
                if self.show_help
                else style(" <yellow>▶</yellow>")
            )
            value_str = a_textwrap(str(value), width=right_col_width)
            lines = value_str.splitlines()
            last_line = lines[len(lines) - 1] if len(lines) else None
            if last_line and len(last_line) + 2 > right_col_width:
                # When the help triangle doesn't fit on the same line.
                show_help_str = f"\n{left_col_filler}" + show_help_str
                help_str_new_line = False
            value_str = style(
                "<editable> "
                + f" </editable>\n{left_col_filler}<editable> ".join(value_str.splitlines())
                + " </editable>"
                + show_help_str,
                nowrap=True,
            )  # noqa
            return value_str

        def __render_helper_text(left_col_filler, right_col_width):
            help_str_lines = style(self.selected_help, nowrap=True).splitlines()
            # When there's line breaks in the help string,
            # we need to apply wrapping per individual line.
            help_str = [a_textwrap(line, width=right_col_width) for line in help_str_lines]
            help_str = "\n".join(help_str)
            first_line_filler = f"\n{left_col_filler} <yellow>" if help_str_new_line else " <yellow>"
            help_str = (
                first_line_filler + f"</yellow>\n{left_col_filler} <yellow>".join(help_str.splitlines()) + "</yellow>"
            )
            help_str = style(help_str, nowrap=True)
            return help_str

        def __render_illegal(illegal_flag, value, left_col_filler, right_col_width):
            illegal_flag_placeholder = "~" * 7  # Just a random character that is unlikely to be repreated in the value.
            value_str = a_textwrap(f"{illegal_flag_placeholder} {value}", width=right_col_width)
            value_str = style(
                "<soft><strikethrough>"
                + f"</strikethrough></soft>\n{left_col_filler} <soft><strikethrough>".join(value_str.splitlines())
                + "</strikethrough></soft>",
                nowrap=True,
            )  # noqa
            value_str = illegal_flag + value_str.replace(illegal_flag_placeholder, "")
            return value_str

        def __render_invalid(value, left_col_filler, right_col_width):
            value_str = a_textwrap(str(value), width=right_col_width)
            value_str = style(
                "<error>" + f"</error>\n{left_col_filler} <error>".join(value_str.splitlines()) + "</error>",
                nowrap=True,
            )  # noqa
            return value_str

        def __render_regular(value, left_col_filler, right_col_width):
            value_str = a_textwrap(str(value), width=right_col_width)
            value_str = f"\n{left_col_filler} ".join(value_str.splitlines())
            return value_str

        def __render_error_msg(left_col_filler, right_col_width):
            error_msg = f"^ {self.inline_errors[self.addresses[index]]}"  # noqa
            error_msg = a_textwrap(error_msg, width=right_col_width)
            error_msg = style(
                left_col_filler
                + "<error_reverse> "
                + f" </error_reverse>\n{left_col_filler}<error_reverse> ".join(error_msg.splitlines())
                + " </error_reverse>",
                nowrap=True,
            )  # noqa
            return error_msg

        # The line_ouput can be wrapping in multiple lines if the value
        # is longer than the available width. We count the lines so we
        # can keep track of the line numbers of our keys.
        def __display_append(line_output):
            for line in line_output.splitlines():
                display.append(line)
                self.line_nr += 1

        #
        #

        _render_one_line(self.data)
        return display

    # Render whatever custom screen is activated (confirm exit, help...).
    def render_active_screen(self):
        return self.active_screen.render()

    ##
    # Editing
    ##

    # Recursively loop through data structure and set value based on key address.
    def value_to_address(self, value, data=None, address=None):
        data = data if data else self.data
        address = address if address else self.selected_address
        reduce(operator.getitem, address[:-1], data)[address[-1]] = value

    # Recursively loop through data structure and remove key based on key address.
    def remove_key_from_address(self, data=None, address=None):
        data = data if data else self.data
        address = address if address else self.selected_address
        reduce(operator.getitem, address[:-1], data).pop(address[-1])

    # Keyboard functionality.
    def process_key_stroke(self, ks):
        self.error = ""

        # For debugging.
        # Log keystroke data.
        # - - - - - - - - - - - - - -
        if 5 in DEBUG:
            self.log(f"\n>> {ks} {ks.name} / {ks.code} / {ks.is_sequence}")
        # - - - - - - - - - - - - - -

        # Schema mode interaction.
        if self.schema_mode:
            if ks.name == "KEY_ESCAPE":
                self.close_schema()
                return

        # Read-only mode.
        if self.read_only:
            if ks.name == "KEY_DOWN":
                self.scroll(1)
                return
            elif ks.name == "KEY_UP":
                self.scroll(-1)
                return

        # Custom screen interaction (confirm exit, help...)
        if self.active_screen and hasattr(self.active_screen, "match_key"):
            key_match = self.active_screen.match_key(ks)
            if key_match:
                return key_match

        # Listen for actions, deactivate listener on invalid input.
        elif self.listen_for_action:
            self.listen_for_action = False
            actions = {
                "?": self.list_commands,
                "p": self.display_path,
                "o": self.open_file,
                "c": self.clear_illegal_keys,
                "f": self.format_file,
                "s": self.open_schema,
                "x": lambda: self.confirm_exit("discard"),
                "KEY_ENTER": lambda: self.confirm_exit("save"),
            }
            if ks in actions:
                self.active_screen = None  # Reset in case command was launched from overview screen.
                exit_msg = actions[ks]()
                return exit_msg
            elif ks.name in actions:
                self.active_screen = None  # Reset in case command was launched from overview screen.
                exit_msg = actions[ks.name]()
                return exit_msg
            else:
                return

        # Edit mode
        elif self.edit_mode:
            # Deactivate action listener.
            self.listen_for_action = False

            # Map key
            # Keys agnostic of edit mode.
            if ks.name == "KEY_UP":
                # Edit mode: move cursor one line up.
                if self.edit_mode:
                    lines = a_textwrap(self.text_input, width=self.right_col_width, drop_whitespace=False).splitlines()
                    lines_len = 0
                    jump_to = 0
                    for i, line in enumerate(lines):
                        is_last_line = i == len(lines) - 1
                        if lines_len + len(line) <= self.cursor_pos and not is_last_line:
                            # Previous lines
                            lines_len += len(line)
                        else:
                            # Current line
                            line_cursor_pos = self.cursor_pos - lines_len
                            prev_line_len = len(lines[max(i - 1, 0)])
                            jump_to = max(lines_len - prev_line_len + min(line_cursor_pos, prev_line_len - 1), 0)
                            break
                    self.cursor_pos = jump_to

                # Default: Go to previous field.
                else:
                    self.store_value()
                    self.go_to_prev_field()
            elif ks.name == "KEY_DOWN":
                # Move cursor one line down.
                lines = a_textwrap(self.text_input, width=self.right_col_width, drop_whitespace=False).splitlines()
                lines_len = 0
                jump_to = 0
                for i, line in enumerate(lines):
                    is_last_line = i == len(lines) - 1
                    if lines_len + len(line) <= self.cursor_pos and not is_last_line:
                        # Previous lines
                        lines_len += len(line)
                    else:
                        # Current line
                        line_cursor_pos = self.cursor_pos - lines_len
                        lines_len += len(line)
                        jump_to = min(lines_len + line_cursor_pos, len(self.text_input))
                        break

                # Move one line down, or jump to end.
                if self.cursor_pos < jump_to:
                    self.cursor_pos = jump_to

                # Go to next field when already at the end.
                else:
                    self.store_value()
                    self.go_to_next_field()
            elif ks.name == "KEY_PGUP":
                # Jump to start.
                self.cursor_pos = 0
            elif ks.name == "KEY_PGDOWN":
                # Jump to end.
                self.cursor_pos = len(self.text_input) - 1
            elif ks.name == "KEY_HOME":
                # Jump to beginning of line.
                lines = a_textwrap(self.text_input, width=self.right_col_width, drop_whitespace=False).splitlines()
                jump_to = 0
                for line in lines:
                    if jump_to + len(line) < self.cursor_pos:
                        jump_to += len(line)
                    else:
                        break
                self.cursor_pos = jump_to
            elif ks.name == "KEY_END":
                # Jump to end of line.
                lines = a_textwrap(self.text_input, width=self.right_col_width, drop_whitespace=False).splitlines()
                jump_to = 0
                for line in lines:
                    # xx += '--' + line
                    if jump_to + len(line) <= self.cursor_pos:
                        jump_to += len(line)
                    else:
                        jump_to += len(line)
                        break
                self.cursor_pos = jump_to - 1
            elif ks.name == "KEY_ENTER":
                self.store_value()
            elif ks.name == "KEY_ESCAPE":
                # Escape key in edit mode is applied with a delay,
                # so we can detect alt+left/right input. This is
                # because alt+left/right comes in as two separate
                # keystrokes.
                def _apply_escape():
                    self.esc_code = False
                    self.cursor_pos = 0
                    self.edit_mode_exit()
                    self.render_one_frame()

                self.esc_code = True
                self.timer = Timer(0.1, _apply_escape)
                self.timer.start()
            elif ks.name == "KEY_BACKSPACE":
                self.text_input = self.text_input[: self.cursor_pos - 1] + self.text_input[self.cursor_pos :]
                self.cursor_pos = max(self.cursor_pos - 1, 0)
            elif ks.name == "KEY_DELETE":
                self.text_input = self.text_input[: self.cursor_pos] + self.text_input[self.cursor_pos + 1 :]
            elif ks.name == "KEY_LEFT":
                if self.cursor_pos > 0:
                    self.set_cursor(-1)
            elif ks.name == "KEY_RIGHT":
                if self.cursor_pos < len(self.text_input):
                    self.set_cursor(1)
            elif str(ks) == "b" and self.esc_code:
                # Alt+left
                self.timer.cancel()
                self.esc_code = False

                # Split words and punctuation.
                words = self.text_input[: self.cursor_pos].split(" ")
                for i, word in list(enumerate(words))[::-1]:
                    if word and word[-1] in ",.;:":
                        words[i] = word[:-1]
                        words.insert(i + 1, word[-1])

                # Jump to start of last words/punctiation.
                last_word = words[len(words) - 1]
                second_last_word = words[len(words) - 2]
                jump_back = len(last_word)
                if jump_back == 0:
                    jump_back += len(second_last_word) + 1
                self.cursor_pos = max(self.cursor_pos - jump_back, 0)
            elif str(ks) == "f" and self.esc_code:
                # Alt+right
                self.timer.cancel()
                self.esc_code = False

                # Split words and punctuation.
                words = self.text_input[self.cursor_pos :].split(" ")
                for i, word in list(enumerate(words))[::-1]:
                    if word and word[-1] in ",.;:":
                        words[i] = word[:-1]
                        words.insert(i + 1, word[-1])

                # Jump to start of next words/punctiation.
                first_word = words[0]
                second_word = words[1] if len(words) > 1 else ""
                jump_fwd = len(first_word)
                if jump_fwd == 0:
                    jump_fwd += len(second_word) + 1
                self.cursor_pos = min(self.cursor_pos + jump_fwd, len(self.text_input))
            elif ks and not ks.is_sequence:
                # Type
                self.type(ks)
            else:
                pass

        # Browse mode
        else:
            if ks.name == "KEY_UP":
                self.go_to_prev_field()
            elif ks.name == "KEY_DOWN":
                self.go_to_next_field()
            elif ks.name == "KEY_RIGHT":
                self.toggle_help(True)
            elif ks.name == "KEY_LEFT":
                self.toggle_help(False)
            elif ks.name == "KEY_ENTER":
                self.edit_mode_enter()
            elif ks.name == "KEY_ESCAPE":
                # Not all screen classes have their own key handlers
                # (eg. all_commands) so for those ESC will always return
                # you to the main screen.
                if self.active_screen:
                    self.active_screen = None
                else:
                    return self.confirm_exit("discard")
            elif ks.name == "KEY_BACKSPACE":
                self.edit_mode_enter(clear=True)
            elif ks == ":":
                self.listen_for_action = True
            elif ks and not ks.is_sequence:
                # Type
                self.edit_mode_enter(clear=True)
                self.type(ks)

    # Enter one character.
    def type(self, ks):
        before_cursor = self.text_input[: self.cursor_pos] + ks
        after_cursor = self.text_input[self.cursor_pos :]
        self.text_input = before_cursor + after_cursor
        self.cursor_pos += 1

    # Enter edit mode.
    def edit_mode_enter(self, clear=False):
        self.edit_mode = True

        # Remove possible field errors.
        if self.selected_address in self.inline_errors:
            self.inline_errors.pop(self.selected_address)

        if not clear:
            # Transfer text input.
            self.text_input = str(self.selected_value)

            # Set cursor at end of line.
            self.set_cursor()

    # Store the text input in self.data.
    def store_value(self):
        data_type = None
        if self.schema:
            data_type_address = self.address_to_schema_address(self.selected_address) + ("type",)
            data_type = self.value_from_address(self.schema, data_type_address)

        if not self.text_input:
            self.text_input = None

        # Format values:
        # - ["a", 'b'] --> ['a', 'b']
        # - [1,2,  3] --> [1, 2, 3]
        # - .1 --> 0.1
        # - etc.
        try:
            self.text_input = ast.literal_eval(self.text_input)
        except BaseException:
            pass

        # Format booleans (this makes boolean input case insensitive).
        if data_type == "boolean" and isinstance(self.text_input, str):
            if self.text_input.lower() == "true":
                self.text_input = True
            elif self.text_input.lower() == "false":
                self.text_input = False

        # Store value in self.data
        self.value_to_address(self.text_input)
        self.edit_mode_exit()

        # Reset
        self.cursor_pos = 0

    # Exit edit mode.
    def edit_mode_exit(self):
        self.validate()
        self.edit_mode = False
        self.text_input = ""

    ##
    # Navigation
    ##

    # Jump to previous field.
    def go_to_prev_field(self):
        prev_index = self.selected_index - 1
        if prev_index < 0:
            prev_index = len(self.addresses) - 1
        self.selected_index = prev_index
        self.selected_line_nr = self.index_to_line_nr_map[self.selected_index]
        self.toggle_help(False)
        self.text_input = ""

    # Jump to next field.
    def go_to_next_field(self):
        next_index = self.selected_index + 1
        if next_index > len(self.addresses) - 1:
            next_index = 0
        self.selected_index = next_index
        self.selected_line_nr = self.index_to_line_nr_map[self.selected_index]
        self.toggle_help(False)
        self.text_input = ""

    # Show/hide helper text (if present in schema).
    def toggle_help(self, show=True):
        self.show_help = show

    # Update cursor.
    # Direction can be 1, -1, or when set to none,
    # the cursor is set to the end of the line.
    def set_cursor(self, direction=None):
        if direction:
            self.cursor_pos += direction
        else:
            # Set cursor at end of line
            self.cursor_pos = len(self.text_input)

    # Scrolling function for read-only mode.
    # Regular mode will scroll automatically
    # via go_to_next_field() & self.selected_line_nr.
    # Search for viewport to find code.
    def scroll(self, jump=1):
        if self.display_height < self.height:
            return
        min = self.height - 1
        max = 0
        self.skip_lines = max(min(self.skip_lines + jump, min), max)  # TODO: Figure out why we need *another* -1 here

    ##
    # Actions
    ##

    # [:?] List all commands.
    def list_commands(self):
        if not self.read_only:
            self.active_screen = self.screens["all_commands"]

    # [:s] Open the schema in the JSON editor.
    def open_schema(self):
        if not self.schema:
            return
        # Store current values
        self.schema_mode_store = {"read_only": self.read_only, "data": self.data}

        self.schema_mode = True
        self.read_only = True
        self.data = self.schema
        self.store_addresses()
        self.set_left_col_width()
        # self.log(self.schema)
        # self.open_json_file(self.schema_path, read_only=True) # %% Trash I think, or useful for linking?

    # Close schema and return to JSON.
    def close_schema(self):
        # for key in self.schema_mode_store:
        #     self[key] = self.schema_mode_store[key] # %% Trash I think, or useful for linking?
        self.schema_mode = False
        self.read_only = self.schema_mode_store["read_only"]
        self.data = self.schema_mode_store["data"]

        self.store_addresses()
        self.set_left_col_width()
        self.validate()

    # Load the previous JSON one. - UNUSED until we implement linked JSON files
    def go_back(self):
        self.options_hist.pop()
        prev_options = self.options_hist[len(self.options_hist) - 1]
        self.open_json_file(*prev_options)

    # [:p] Display the full JSON path.
    def display_path(self):
        self.show_path = not self.show_path

    # Check if JSON matches the schema.
    # In strict mode, all missing fields will be counted,
    # otherwise only required fields will be counted.
    def detect_discrepancies(self):
        def _recursion(data_item, schema_props, required_keys):
            nonlocal missing_keys
            for key in schema_props:
                # Detect if it's a parent field.
                is_dict = False
                if "properties" in schema_props[key]:
                    is_dict = True

                # Check if key is required.
                is_required = key in required_keys

                # Add missing values.
                if key not in data_item and (self.strict or is_required):
                    self.log(is_required, key)
                    return True

                # Continue recursively.
                if is_dict and key in data_item:
                    required_keys = schema_props[key]["required"] if "required" in schema_props[key] else []
                    return _recursion(data_item[key], schema_props[key]["properties"], required_keys)

        #
        #

        required_keys = self.schema["required"] if (self.schema and "required" in self.schema) else []
        missing_keys = (
            _recursion(self.data, self.schema["properties"], required_keys)
            if (self.schema and "properties" in self.schema)
            else False
        )
        illegal_keys = len(self.illegal_errors)

        return missing_keys, illegal_keys

    # [:f] Reformat JSON to match schema.
    def format_file(self):
        if not self.schema:
            return

        # Format one level.
        def _recursion(data, schema_props):
            for key in schema_props:
                # Detect if it's a parent field.
                is_dict = False
                is_array = "type" in schema_props[key] and schema_props[key]["type"] == "array"
                if not is_array:
                    if "properties" in schema_props[key]:
                        is_dict = True

                # Add missing values.
                if key not in data:
                    val = {} if is_dict else schema_props[key]["default"] if "default" in schema_props[key] else ""
                    data[key] = val

                    # Maybe we only want to reset required fields?
                    # is_required = 'required' in data[key] and data[key]['required']
                    # if is_required:
                    #     data[key] = val

                # Continue recursively.
                if is_dict:
                    _recursion(data[key], schema_props[key]["properties"])

        #
        #

        self.clear_illegal_keys()
        _recursion(self.data, self.schema["properties"])
        self.data = self.reorder_by_schema()
        self.store_addresses()
        self.set_left_col_width()
        self.validate()

    # Reorder data accoridng to schema order.
    # Note: schema data needs to match json data.
    def reorder_by_schema(self):
        ordered_data = OrderedDict()

        def _recursion(data, schema, ordered_data):
            if "properties" in schema:
                for key, val in schema["properties"].items():
                    if isinstance(val, dict) and "properties" in val:
                        ordered_data[key] = {}
                        _recursion(data[key], val, ordered_data[key])
                    else:
                        ordered_data[key] = data[key]

        _recursion(self.data, self.schema, ordered_data)
        return ordered_data

    # [:c] Remove illegal fields from self.data.
    def clear_illegal_keys(self):
        for address in self.illegal_errors:
            self.remove_key_from_address(address=address)
        self.illegal_errors = {}

    # [:o] Open file in text editor.
    def open_file(self):
        import subprocess

        subprocess.run(["open", self.json_path])

    # [:ENTER] [ESC] [ctrl+c] Confirm discarding of changes & exit.
    # - - -
    # [ctrl+c] --> [DISCARD] preselected
    # [ESC] --> [CANCEL] preselected
    # [:ENTER] --> [SAVE] preselected
    def confirm_exit(self, sel="save"):
        # Check if there has been any changes.
        def _has_changes():
            try:
                with open(self.json_path, "r", encoding="UTF-8") as f:
                    file_data = json.load(f)
                if self.data == file_data:
                    return False
                else:
                    return True
            except BaseException:
                return True

        #
        #

        if _has_changes():
            self.active_screen = self.screens["confirm_exit"]
            self.active_screen.sel = sel
        elif sel == "cancel":
            # When prssing ESC without making any changes.
            return style("<soft>No changes were made.</soft>")
        else:
            return True

    ##
    # Utility
    ##

    # Recursively loop through data structure and select
    # value based on an address (tuple with key tree).
    def value_from_address(self, data, address):
        try:
            result = reduce(lambda d, k: d.get(k) if isinstance(d, dict) else d[int(k)], address, data)
        except BaseException as err:
            result = None

        # Do another pass in case the schema is referenced.
        # This is basically doing the same as the reduce but
        # this time we expand referenced schemas (lookup_schema_def).
        if data == self.schema and not result:

            def _recursion(schema_prop, addy):
                key = addy.popleft()

                # Abort for illegal values.
                if key not in schema_prop:
                    return

                val = schema_prop[key]
                if isinstance(val, dict) and len(addy):
                    return _recursion(val, addy)
                else:
                    return val

            #
            #

            result = _recursion(data, deque(address))

        return result

    # Convert a json address into the equivalent schema address.
    def address_to_schema_address(self, address):
        schema_address = []
        for item in address:
            schema_address.extend(("properties", item))
        return tuple(schema_address)

    # Returns the total length of the ANSI escape codes.
    # Useful when we want to stuff an ANSI code containing
    # string with filler characters (:<) and we need to know
    # by how much we need to extend the filler length.

    def len_ansi(self, text):
        pattern_ansi = r"\x1b\[[0-9;]*[m]"
        ansi_codes = re.findall(pattern_ansi, text)
        total_length = sum(len(code) for code in ansi_codes)
        return total_length

    ##
    # Development
    ##

    # Print any debug data below the JSON content.
    # Usage: self.log()
    def log(self, *args, override=False):
        args_str = [str(arg) for arg in args]
        args_str = " / ".join(args_str)
        if override:
            self.log_values = [args_str]
        else:
            self.log_values.append(args_str)


# exit_screen_show
class ConfirmExitScreen:
    name = "confirm_exit"
    keep_title = True
    actions = ["save", "discard", "cancel"]
    sel = None  # Which exit button gets pre-selected from self.exit_screen_actions

    def __init__(self, edit_json_class):
        self.ej = edit_json_class

    def render(self):
        actions_off = [f" {action.upper()} " for action in self.actions]
        actions_on = [f"<reverse> {action.upper()} </reverse>" for action in self.actions]
        output = []
        for i, action in enumerate(self.actions):
            if action == self.sel:
                output.append(actions_on[i])
            else:
                output.append(actions_off[i])
        output = " ".join(output)
        return style(output + "  ◀ ▶", nowrap=True)

    def match_key(self, ks):
        if ks.name == "KEY_LEFT":
            index = self.actions.index(self.sel)
            index = max(index - 1, 0)
            self.sel = self.actions[index]
        elif ks.name == "KEY_RIGHT":
            index = self.actions.index(self.sel)
            index = min(index + 1, len(self.actions) - 1)
            self.sel = self.actions[index]
        elif ks.name == "KEY_ESCAPE":
            return self.exit()
        elif ks.name == "KEY_ENTER":
            return self.submit()
        elif ks == "y":
            return "<yellow>Your changes have been discarded!!</yellow>"
        else:
            return

    # Actions

    # Cancel & return to editor.
    def exit(self):
        self.ej.active_screen = None

    # Submit exit screen choice.
    def submit(self):
        actions = {
            "save": self.save_changes,
            "discard": self.discard_changes,
            "cancel": self.exit,
        }
        action = actions[self.sel]
        exit_message = action()
        return exit_message

    # Exit & save changes.
    def save_changes(self):
        self.exit()
        if bool(self.ej.inline_errors):
            self.ej.error = "Some of your changes are invalid."
        else:
            try:
                # raise ValueError("Fake error for debugging")
                if self.ej.demo:
                    return style(
                        "<green>Changes saved to JSON file.</green>\n<soft>Note: you were in demo mode so your changes have not really been saved.</soft>",
                        pad=1,
                    )
                    # return self.ej.term.green('Changes saved to JSON file.') + self.ej.term.yellow('\nNote: you were in demo mode so your changes have not really been saved.')  # noqa
                else:
                    with open(self.ej.json_path, "w", encoding="UTF-8") as f:
                        json.dump(self.ej.data, f, indent=4)
                    # return self.ej.term.green('Changes saved to JSON file.')
                    return style("<green>Changes saved to JSON file.</green>", pad=1)
            except ValueError as err:
                self.ej.error = "There was an error saving your changes."

    # Exit & discard changes.
    def discard_changes(self):
        self.exit()
        # self.exit_screen_show = False
        return "<yellow>Your changes have been discarded.</yellow>"


class AllCommandsScreen:
    name = "all_commands"
    keep_title = False
    commands = {
        "?": "List all commands",
        "p": "Show/hide JSON file path",
        "o": "Open JSON file in text editor",
        "c": "Clear any illegal fields",
        "f": "Format JSON file to match schema",
        "s": "Open schema file",
        "x": "Discard changes",
        "ENTER": "Save changes",
    }

    def __init__(self, edit_json_class):
        self.ej = edit_json_class

    def render(self):
        output = []
        output.append(style(f"<yellow>Available Commands</yellow>", nowrap=True))
        output.append(self.ej.sep)
        output.append(style(f"<reverse> ESC </reverse> Go back", nowrap=True))
        output.append(self.ej.sep)
        for key, action in self.commands.items():
            if key == "c" and (not self.ej.illegal_errors or not self.ej.schema):
                action = f"<soft>{action}</soft>"
            if (key == "f" or key == "s") and not self.ej.schema:
                action = f"<soft>{action}</soft>"
            output.append(style(f"<cmd>:{key:<{7}}</cmd> {action}", nowrap=True))
        return "\n".join(output)


if __name__ == "__main__":
    ej = EditJson()
    ej(demo=True)
