# Creating Your Own Toolkit <!-- omit from toc -->

Integrating your own workflows into OpenAD is relatively straightforward.

### Table of Contents <!-- omit from toc -->
- [Setup](#setup)
  - [metadata.json](#metadatajson)
  - [login.py](#loginpy)
- [Adding functions](#adding-functions)
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

Each command is contained within a single JSON file that follows the structure "func_<funcname>.json". "func_hello_world.json".

The structure looks as follows:

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
        The category under which this command is organized. If you don't have different categories for yout commands, you can set this to "General".
    -   `command`
        The structure of the command as it will be displayed in the list of commands, or when a user requests help about this particular command. See [Command Documentation](#command-documentation) below for guidance about notation.
    -   `description`<br>
        A description of what the command does. If your command description is going to be longer than one or two lines, we recommand to set this value to an empty string "" and store the description in a separate text file with the same base filename with "--decription.txt" appended. So if your function file were to be called "func_hello_world.json", your description file would be called "func_hello_world--description.txt".

        > **Important:** Because the desciption text is parsed for different outputs (CLI, Jupyter, HTML, API) it needs to use a specific OpenAD styling syntax, the documentation for which you can find [here].
    -   `url`<br>
        TBD - A link to online documentation for this command.
-   `library`<br>
    TBD
-   `method`<br>
    TBD
-   Command clauses: TBD
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

### Command Documentation


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

<br>

### Publishing your toolkit

Once your toolkit adheres to the aformentioned specifications, [get in touch](https://acceleratedscience.github.io/openad-docs/about.html) so we can review it and add it to the publicly available OpenAD tookits.