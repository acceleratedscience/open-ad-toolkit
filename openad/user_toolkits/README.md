> **WARNING:** This documentation page is still under construction. It is incomplete and may have incorrect information.

---

# Creating Your Own Toolkit <!-- omit from toc -->

Integrating your own workflows into OpenAD is relatively straightforward.

### Table of Contents <!-- omit from toc -->
- [Setup](#setup)
  - [metadata.json](#metadatajson)
  - [login.py](#loginpy)
- [Adding Commands](#adding-commands)
  - [1. func\_hello\_world.json](#1-func_hello_worldjson)
  - [2. func\_hello\_world.py](#2-func_hello_worldpy)
  - [3. func\_hello\_world.txt](#3-func_hello_worldtxt)
  - [Command Notation](#command-notation)
- [Publishing a Toolkit](#publishing-a-toolkit)
  - [description.txt](#descriptiontxt)
  - [oneline\_desc.txt](#oneline_desctxt)
  - [Publishing your toolkit](#publishing-your-toolkit)

<br><br>

## Setup

The toolkit architecture depends on a few basic files to work. You can copy the [DEMO](./DEMO) toolkit directory to hit the ground running.


1. [`metadata.json`](#metadatajson)
2. [`login.py`](#loginpy) (optional)

<br>

### metadata.json

This mandatory file is responsible for generating the toolkit's splash screen. This is displayed when the user enters the toolkit context by running `set context <toolkit_name>`.

    {
        "banner": "DEMO",
        "title": "This is a Demo Toolkit",
        "author": "Jane Doe",
        "version": "0.0.1",
        "intro": "This toolkit is meant as a demonstration on how to set up a toolkit. This intro paragraph should contain a brief description of what the toolkit does, ideally not much longer than ~250 characters.",
        "commands": {
            "hello world": "Say hello",
            "demo docs": "Read the docs in your browser"
        }
    }

The splash screen generated from the JSON file above looks like this:

![demo-splash-page](readme/demo-splash.png)

<br>

### login.py

Here you can expose your authentication API. If this file is present, is will be called whenever the user enters the toolkit context by running `set context <toolit_name>`. If this file is not present, authentication will be skipped.
<!-- Todo: remove need for login file -->

The [`login.py`](./DEMO/login.py) template takes care of success and error handling and ensures a unified user experience across all tookits. Instructions are in the file. You may have to customize it a bit more if your authentication API doesn't follow jwt/host/email/api_key conventions.

<br><br>

## Adding Commands

Each command is contained in two (or three) files. We'll document the "hello world" example from our demo toolkit.
1. [`func_hello_world.json`](#func_hello_worldjson) The command file (JSON) which contains the command structure and help parameters.
2. [`func_hello_world.py`](#func_hello_worldpy) The function file (Python) which contains the command function to be executed.
3. Optional: [`func_hello_world.txt`] The description file (text) which may contain a more elaborate command description, when a field in the JSON command file is not enough.

Per function, these three files should be stored in the same directory. However the total of all your functions can be organized into a directory structure as desired. We will scan your entire toolkit directory tree and parse any command files we find. 

Multiple JSON command files may point to a single Python function file containing multiple methods. They are linked through the `library` parameter in the JSON command file.

On the other hand, only one description text file can exist per JSON command file, and they are linked by having the same base filename.

<br>

### 1. func_hello_world.json

The command file's filename follows the structure `func_\<func_name\>.json`. Hence our "hello world" command file would be named `func_hello_world.json`.

Our `func_hello_world.json` file structure would look as follows:

    {
        "fplugin": "demo",
        "command": "hello world", 
        "subject": "hello world",
        "exec_type": "standard_statement",
        "exec_cmd": "hello world",
        "help": {
            "name": "hello world",
            "category": "General",
            "command": "hello world",
            "description": "Display \"Hello, world\"",
            "url": "https://helloworld.app/docs#hello-world"
        },
        "library": "hello_world",
        "method": "print_hello_world",
        (Command clauses)
    }

<details>
<summary><b>JSON Breakdown</b></summary>

-   `fplugin`<br>
    The name of your toolkit, the same as your tookit directory (case insensitive).
-   `command`<br>
    TBD
-   `subject`<br>
    TBD
-   `exec_type`<br>
    TBD "standard_statement" or "method"
-   `exec_cmd`<br>
    TBD
-   `help`<br>
    Contains all we need to plug into the OpenAD help functionality.
    -   `name`<br>
        TBD
    -   `category`<br>
        The category under which this command is organized. If you don't have different categories for your commands, you can set this to "General".
    -   `command`
        The structure of the command as it will be displayed in the list of commands, or when a user requests help about this particular command. See [Command Documentation](#command-documentation) below for guidance about notation.
    -   `description`<br>
        A description of what the command does. If your command description is going to be longer than one or two lines, we recommand to set this value to an empty string "" and store the description in a separate text file with the same base filename with "--decription.txt" appended. So if your function file were to be called `func_hello_world.json`, your description file would be called `func_hello_world--description.txt`.

        > **Important:** Because the desciption text is parsed for different outputs (CLI, Jupyter, HTML, API) it needs to use a specific OpenAD styling syntax, the documentation for which you can find [here].
    -   `url`<br>
        TBD - A link to online documentation for this command.
-   `library`<br>
    The base name (no extension) of the Python file in the same directory as the JSON file, that contains the command function.
-   `method`<br>
    TBD

</details>

<br>

<details>
<summary><b>Command Clauses</b></summary>

These are common built-in command patterns that represent certain behaviors. Together with the `command` <!-- @Phil -->parameter of your JSON file, they define the structure of your command.
-   `SAVE_AS`<br>
    This (always optional) clause is meant for functions that output data, and should cause the output of your command function to be saved to disk instead of being displayed.
    
    -   **JSON notation:** `"SAVE_AS": {}`
    -   **Function access:** `if "save_as" in inputs:`
    -   **Command notation:** `hello world [ save as '<filename.txt>' ]`

-   `ESTIMATE_ONLY`<br>
    This (always optional) clause is meant for functions that may take a long time to return results, and should cause your function to return an estimate of result count rather than the actual results.

    -   **JSON notation:** `"ESTIMATE_ONLY": {}`
    -   **Function access:** `if "estimate_only" in inputs:`
    -   **Command notation:** `hello world [ estimate only ]`
-   `RETURN_AS_DATA`<br>
    This (always optional) clause is meant for functions that return styled data, and should remove any styling from your data so it can be consumed by endpoints where the styling is not welcome.

    -   **JSON notation:** `"RETURN_AS_DATA": {}`
    -   **Function access:** `if "return_as_data" in inputs:`
    -   **Command notation:** `hello world [ return as data ]`
-   `SINGLE_PARM`<br>
    TBD
    
    -   **JSON notation:** `"SINGLE_PARM": { "smiles": "desc" }`
    -   **Function access:** `if "return_as_data" in inputs:`
    -   **Command notation:** `hello world [ return as data ]`
-   `SHOW`<br>
    This (always optional) clause is meant for functions that may return different types or formats of data, and should be used to specify what kind of data to return.

    -   **JSON notation:** `"SHOW": ["foo", "bar"]`
    -   **Function access:**
        ```
        if "show_data" in inputs:
            for option in inputs["show_data"]:
                if option == "foo":
        ```     
    -   **Command notation:** `hello world show (foo | bar)`
-   `USING`<br>
    This clause is meant to pass a number of custom, optional variables to your function.

    -   **JSON notation:**
        ```
        "USING": {
            "id": "str",
            "page_size": "int",
            "foo_bar": "desc"
        }
        ```
        <!-- @Phil desc? -->
    -   **Function access:**
        ```
        if "page_size" in inputs:
            page_size = inputs["page_size"][val]
        ```
    -   **Command notation:** `hello world using (id='<str>' page_size=<int> foo_bar=?)`
-   `FROM`<br><!-- @Phil needs to be made uppercase -->
     This clause is meant for functions that allow for different input types, where it indicates the provided input type. Only "dataframe", "file" and "list" are supported as input types.

    -   **JSON notation:** `"FROM": ["dataframe", "file", "list"]`
    -   **Function access:**
        ```
        if isinstance(inputs["from_source"], dict) and inputs["from_source"]["from_list"] is not None:
            from_list = inputs["from_source"]["from_list"]
        elif "from_list" in inputs["from_source"][0]:
            from_list = inputs["from_source"][0]["from_list"]
        elif "from_dataframe" in inputs:
            dataframe = inputs["from_dataframe"]
        elif "from_file" in inputs:
            from_file = inputs["from_file"]
        ```
        <!-- @Phil why not "from_list" in inputs? -->
    -   **Command notation:** `hello world from dataframe <dataframe_name> | file '<csv_filename>' | list ['<string>','<string>']`
-   `USE_SAVED`<br><!-- @Phil needs to be made uppercase -->
    This clause is meant for functions that use cacheing, where it indicated that a cashed result can be used instead of re-running the function.

    -   **JSON notation:** `"USE_SAVED": "True"` <!-- @Phil why not {}? func_predict_retro has "use_saved": "False" -->
    -   **Function access:** `if "use_saved" in inputs:`
    -   **Command notation:** `hello world [ use_saved ]`
-   `ASYNC`<br><!-- @Phil needs to be made uppercase -->
    TBD

    -   **JSON notation:** `"ASYNC": "both"` `"ASYNC": "only"`
    -   **Function access:** `TBD (do_async/async)`
    -   **Command notation:** `TBD`

<!--
- "SAVE_AS": {}
- "SINGLE_PARM": { "collection": "desc" }
- "SINGLE_PARM": { "search_string": "desc" },
- "SINGLE_PARM": { "domain": "desc" }
"USING": {
    "page_size": "int",
    "system_id": "str",
    "edit_distance": "int",
    "display_first": "int"
},
- "SHOW": ["DATA", "DOCS"],
- "ESTIMATE_ONLY": {},
- "SAVE_AS": {},
- "RETURN_AS_DATA": {}
- "SINGLE_PARM": { "smiles": "desc" },
x "OLD_USING": { "prediction_id": "desc", "ai_model": "desc", "topn": "int" },
x "USING": { "ai_model": "desc", "topn": "int" },
x "USING_old": { "prediction_id": "desc", "ai_model": "desc" },
x "USING": { "prediction_id": "desc", "ai_model": "desc" },
- "SINGLE_PARM": { "molecule": "desc" },
x "USING": { "prediction_id": "desc", "ai_model": "desc" },
- "SINGLE_PARM": { "molecule": "desc" },
x "USING": {
    "availability_pricing_threshold": "int",
    "available_smiles": "desc",
    "exclude_smiles": "desc",
    "exclude_substructures": "desc",
    "exclude_target_molecule": "str",
    "fap": "number",
    "max_steps": "int",
    "nbeams": "int",
    "pruning_steps": "int",
    "ai_model": "desc"
},
- "SINGLE_PARM": { "recipe": "desc" },
-->

</details>

<br>

### 2. func_hello_world.py
<!-- @Phil thoughts on keeping the filename consistent? -->

This file contains your command function. It gets passed two parameters:
- `inputs`: A dictionary containing all information relating to the command the user typed.
- `cmd_pointer`: An instance of the [RUNCMD class](../app/main.py), which is a subclass of the [`Cmd`](https://docs.python.org/3/library/cmd.html) library class.

What our "hello world" example's function file looks like:

    def hello_world(inputs: dict, cmd_pointer):
        print('hello, world')
        # to do: print inputs and cmd_pointer content

What it outputs:

    hello, world

    inputs:
    - ...

    cmd_pointer:
    - ...

<br>

### 3. func_hello_world.txt

Lorem

<br>

### Command Notation

-   Optional clauses should be encapsulated in square brackets padded with a space.

        hello world [ repeat 2 ]
    
-   Variable names should be displayed with angle brackets and underscores instead of spaces. When a variable is to be quoted, make sure to include the quotation marks in the command.

        hello [ <audience_name> ]
        hello [ '<first_name>' ]

-   When describing different options for a clause, list them separated by a pipe. For long commands, it may be unclear which word is part of the main command or the clause options. To avoid confusion, make sure to add examples to your command documentation that will resolve any confusion.

        hello world pink | red | green

-   When describing a filename, add the extension(s) in the variable name, as such:

        hello world [ save as <filename.txt> ]
        hello world [ save as <csv_or_sdf_filename> ]

- When describing lists, make sure to notate them without spaces between the square brackets, to avoid confusion with optional clauses. Also make it clear what is supposed to go in the list. Use ellipsis to indicate whether the length of the list can be infinite.

        hello ['<first_name>','<first_name>',...]



<br><br>

## Publishing a Toolkit

By publishing a toolkit to the OpenAD comminity, it will be made available for all OpenAD users to install. This means that it will be displayed in the list of toolkits when you run `list all toolkits`, and users will be able to install it simply by running `add toolkit <toolkit_name>`.

If you wish to publish your toolkit, a few more files are required.
1. [`description.txt`](descriptiontxt)
2. [`oneline_desc.txt`](oneline_desctxt)

<br>

### description.txt

The `description.txt` file is used to train the LLM with the toolkit functionality. It should contain a detailed, unambiguous description of how your toolkit works and what it is meant to achieve. You can look at the other toolkits for inspiration.

At the bottom of your file, on a separate line, you should include the following line, verbatim:

    The following commands are available for this toolkit:

Then you should run the script below, which gathers all your toolkit commands and lists them at the bottom of the description file. The script will look for the line described above and replace everything after with the updated commands. If this exact line is not present, the script will abort and throw an error.

    python openad/user_toolkits/<toolkit_name>/description_update.py

<br>

### oneline_desc.txt

This file contains a very brief description of the toolkit, using only 4-5 words. This wil be displayed when listing available toolkits.

![toolkit-list](readme/toolkit-list.png)

<br>

### Publishing your toolkit

Once your toolkit adheres to the aformentioned specifications, [get in touch](https://acceleratedscience.github.io/openad-docs/about.html) so we can review it and add it to the publicly available OpenAD tookits.