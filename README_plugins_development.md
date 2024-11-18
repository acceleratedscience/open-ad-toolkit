<sub>[&larr; BACK](./README.md#openad)</sub>

# OpenAD Plugin Development Manual

Creating a plugin for OpenAD is simple. Here's how you do it.

1. **Clone the `plugin_kit` repository**

        git clone git@github.com:acceleratedscience/openad-plugin-kit.git

1. **Provide the Plugin Parameters**
    
    Rename the plugin folder as desired, let's say `my_plugin`.

    Locate `my_plugin/plugin_params.py` and set the required parameters:

    - `LIBRARY_NAME`: The name of your plugin as it is displayed to OpenAD users. It is also how users can list all your plugin's commands: `plugin demo ?`
    - `LIBRARY_KEY`: A snake_case version of your plugin's name, used to create internal command ids.
    - `PLUGIN_NAMESPACE`: A short word or letter combination that will prepend all your commands.
    - `CMD_NOTE`: An optional information note that can be added to all your commands, for example: `To see all available demo comands, run 'demo plugin ?'`

1. **Install the plugin in edit mode**
    
    This lets you test your plugin in OpenAD while editing it.

        cd /path/to/my_plugin
        pip install . -e
    
    Reboot OpenAD for any changes to take effect.

        reboot

1. **Copy the hello_world command folder**

    Rename it to a snake_case name for your command, let's say `my_first_command`.

1. **Find the plugin file: `/my_plugin/my_first_command/command.py`**

    Rename this file to a snake_case name for your command and insert your own functionality. The file is self-explanatory.

1. **Build your command**
   
   The command is constructed by `statements.append()` using [pyparsing](https://github.com/pyparsing/pyparsing/).
   
   Every word in your command will be represented by a pyparsing object in `/my_plugin/command/plugin_grammar.py`.
   
   For regular words, this is as simple as `py.CaselessKeyword('hello')` while single characters are represented as `py.Literal('*')`. For more advanced types like lists, molecule identifiers, optional quotes and more, we provide a ready-made library with pyparsing building blocks:

   [GRAMMAR DEFINITIONS](https://github.com/acceleratedscience/open-ad-toolkit/tree/main/openad/core/grammar_def.py) 

   For more advanced requirements, you can find the list of pyparsing types in the [pyparsing documentation](https://pyparsing-docs.readthedocs.io/en/latest/pyparsing.html).

1. **Publish your repository on GitHub**

1. **Install your plugin**
    
    First uninstall the developemnt version.

        pip uninstall my_plugin
    
    Then install the plugin from GitHub.

        pip install git+https://github.com/my_username/my_plugin.git