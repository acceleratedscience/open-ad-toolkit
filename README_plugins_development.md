<sub>[&larr; BACK](./README.md#openad)</sub>

# OpenAD Plugin Development Manual

A simple manual on how to create your first OpenAD plugin in minutes.

A plugin is an easy way for your application to integrate with all of OpenAD's functionality like the molecule viewer and more.

<br>

### 1. Clone the _plugin_kit_ repository

    git clone git@github.com:acceleratedscience/openad-plugin-kit.git
    
```plaintext
my_plugin/
├── main.py
├── plugin_grammar.py
├── plugin_params.py
├── commands/
│   ├── hello_world/
│   │   ├── command.py
│   │   ├── description.txt
│   │── ...
├── tests/
│   ├── test_hello_world.py
│   ├── ...
└── README.md
```


<br>

### 2. Provide the Plugin Parameters
    
Rename the plugin folder as desired. Currently it's called `my_plugin`.

Locate `my_plugin/plugin_params.py` and set the required parameters:

- **`LIBRARY_NAME`**<br>
    The name of your plugin as it is displayed to OpenAD users.<br>
    It is also how users can list all your plugin's commands: `plugin demo ?`<br>
    <br>
- **`LIBRARY_KEY`**<br>
    A snake_case version of your plugin's name, used to create an internal identifier.<br>
    Avoid using the word "plugin" here to prevent redundancy.<br>
    <br>
- **`PLUGIN_NAMESPACE`**<br>
    A short word or letter combination that will prepend all your commands.<br>
    <br>
- **`CMD_NOTE`**<br>
    An optional informational note that can be added to all your commands, for example:<br>
    `To see all available demo comands, run 'plugin demo ?'`

<br>

### 3. Install the plugin in edit mode
    
This lets you test your plugin in OpenAD while editing it.

    cd /path/to/my_plugin
    pip install -e .

Reboot OpenAD for any changes to take effect.

    reboot

<br>

### 4. Copy the hello_world command folder

Rename it to a snake_case name for your command, let's say `my_first_command`.

<br>

### 5. Find the plugin file: `/my_plugin/my_first_command/command.py`

Insert your own functionality. The file is self-explanatory.

<br>

### 6. Build your command
   
The command is constructed by `statements.append()` using [pyparsing](https://github.com/pyparsing/pyparsing/).

Every word in your command will be represented by a pyparsing object in `/my_plugin/plugin_grammar.py`.

For regular words, this is as simple as `py.CaselessKeyword('hello')` while single characters are represented as `py.Literal('*')`. For more advanced types like lists, molecule identifiers, optional quotes and more, we provide a ready-made library with pyparsing building blocks:

[GRAMMAR DEFINITIONS](https://github.com/acceleratedscience/open-ad-toolkit/tree/main/openad/core/grammar_def.py) 

For more advanced requirements, you can find the list of pyparsing types in the [pyparsing documentation](https://pyparsing-docs.readthedocs.io/en/latest/pyparsing.html).

<br>

### 7. Edit `pyproject.toml` to reflect your plugin

You may also want to update the `README.md`

<br>

### 8. Publish your repository on GitHub

<br>

### 9. Install your plugin
    
First uninstall the developemnt version.

    pip uninstall my_plugin

Then install the plugin from GitHub.

    pip install git+https://github.com/my_username/my_plugin.git