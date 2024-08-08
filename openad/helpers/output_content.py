"""
The introductionary text generated for the `intro` command.
"""

from openad.helpers.files import open_file

# Read description files
description, err_msg = open_file("docs/source/description.txt", return_err=True)
if not description:
    description = "<error>(openAD description not found)</error>"
about_workspace, err_msg = open_file("docs/source/about_workspace.txt", return_err=True)
if not about_workspace:
    about_workspace = "<error>(Workspace description not found)</error>"
about_plugin, err_msg = open_file("docs/source/about_plugin.txt", return_err=True)
if not about_plugin:
    about_plugin = "<error>(Plugin description not found)</error>"
about_context, err_msg = open_file("docs/source/about_context.txt", return_err=True)
if not about_context:
    about_context = "<error>(Context description not found)</error>"
about_run, err_msg = open_file("docs/source/about_run.txt", return_err=True)
if not about_run:
    about_run = "<error>(Run description not found)</error>"


openad_intro: str = f"""<h1>Introduction to OpenAD</h1>

{description}

<link>accelerate.science/projects/openad</link>



<h2>Workspaces</h2>
{about_workspace}

<soft>To see how to work with workspaces, run <cmd>? workspace</cmd></soft>



<h2>Plugins</h2>
{about_plugin}

<soft>To see how to work with plugins, run <cmd>? toolkit</cmd></soft>



<h2>Context</h2>
{about_context}

<soft>To see how to change the context, run <cmd>? context</cmd></soft>



<h2>Runs</h2>
{about_run}

<soft>To see how to work with runs, run <cmd>? run</cmd></soft>"""
