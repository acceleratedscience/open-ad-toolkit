# Creating Your Own Toolkit <!-- omit from toc -->

Integrating your own workflows into OpenAD is relatively straightforward.

### Table of Contents <!-- omit from toc -->
- [Setup](#setup)
  - [metadata.json](#metadatajson)
  - [login.py](#loginpy)
- [Adding functions](#adding-functions)
  - [Command clauses](#command-clauses)
  - [Command Documentation](#command-documentation)
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

Here you can expose your authentication API. If this file is present, is will be called whenever the user enters the toolkit context by running `set context <toolit_name>`. If this file is not present, authentication will be skipped. (@Phil correct?)

The [`login.py`](./DEMO/login.py) template takes care of success and error handling and ensures a unified user experience across all tookits. Instructions are in the file. You may have to customize it a bit more if your authentication API doesn't follow jwt/host/email/api_key conventions. (@Phil what about other variables like username instead of email?)

<br><br>

## Adding functions

Each command is contained within a single JSON file whose name follows the structure `func_\<funcname\>.json`. A hello world function file would thus be called `func_hello_world.json`.

You can store your function files flat, or organize them into a directory structure as desired. We will scan your entire toolkit directory tree for function files and parse them all.

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
    TBD
-   `method`<br>
    TBD

<br>

### Command clauses

These are common built-in command patterns that represent certain behaviors. Together with the `command` (@ph) parameter of your JSON file, they define the structure of your command.
-   `SAVE_AS`<br>
    This (always optional) clause is meant for functions that output data, and should cause the output of your command function to be saved to disk instead of being displayed.
    
    - **JSON notation:** `"SAVE_AS": {}`
    - **Function access:** `if "save_as" in inputs:`
    - **Command notation:** `hello world [ save as '<filename.txt>' ]`

-   `ESTIMATE_ONLY`<br>
    This (always optional) clause is meant for functions that may take a long time to return results, and should cause your function to return an estimate of result count rather than the actual results.

    - **JSON notation:** `"ESTIMATE_ONLY": {}`
    - **Function access:** `if "estimate_only" in inputs:`
    - **Command notation:** `hello world [ estimate only ]`
-   `RETURN_AS_DATA`<br>
    This (always optional) clause is meant for fuctions that return styled data, and should remove any styling from your data so it can be consumed by endpoints where the styling is not welcome.

    - **JSON notation:** `"RETURN_AS_DATA": {}`
    - **Function access:** `if "return_as_data" in inputs:`
    - **Command notation:** `hello world [ return as data ]`
-   `SINGLE_PARM`<br>
    TBD
    
    - **JSON notation:** `"SINGLE_PARM": { "smiles": "desc" }`
    - **Function access:** `if "return_as_data" in inputs:`
    - **Command notation:** `hello world [ return as data ]`
-   `SHOW`<br>
    This (always optional) clause is meant for fuctions that may return different types or formats of data, and should be used to specify what kind of data to return.

    - **JSON notation:** `"SHOW": ["foo", "bar"],`
    - **Function access:**
        
        if "show_data" in inputs:
            for option in inputs["show_data"]:
                if option == "foo":
                
    - **Command notation:** `hello world show(foo | bar)`
-   `AAA`<br>
    Foo
-   `AAA`<br>
    Foo
-   `AAA`<br>
    Foo
-   `AAA`<br>
    Foo
-   `AAA`<br>
    Foo
-   `AAA`<br>
    Foo
-   `AAA`<br>
    Foo

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
"SHOW": ["DATA", "DOCS"],
- "ESTIMATE_ONLY": {},
- "SAVE_AS": {},
- "RETURN_AS_DATA": {}
- "SINGLE_PARM": { "smiles": "desc" },
x "OLD_USING": { "prediction_id": "desc", "ai_model": "desc", "topn": "int" },
"USING": { "ai_model": "desc", "topn": "int" },
x "USING_old": { "prediction_id": "desc", "ai_model": "desc" },
x "USING": { "prediction_id": "desc", "ai_model": "desc" },
- "SINGLE_PARM": { "molecule": "desc" },
x "USING": { "prediction_id": "desc", "ai_model": "desc" },
- "SINGLE_PARM": { "molecule": "desc" },
"USING": {
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

<br>

### Command Documentation

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

@Phil how can do `from docs.generate_docs import render_description_txt`

<br>

### oneline_desc.txt

This file contains a very brief description of the toolkit, using only 4-5 words. This wil be displayed when listing available toolkits.

![toolkit-list](readme/toolkit-list.png)

<br>

### Publishing your toolkit

Once your toolkit adheres to the aformentioned specifications, [get in touch](https://acceleratedscience.github.io/openad-docs/about.html) so we can review it and add it to the publicly available OpenAD tookits.