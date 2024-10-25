"""
The introductionary text generated for the `intro` command.
"""

import os
import sys
from openad.helpers.files import open_file

# Add the root directory to the sys.path
root_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
if str(root_dir) not in sys.path:
    sys.path.append(root_dir)
# for path in sys.path:
#     print("#", path)

# Construct the file paths relative to the script directory
description_path = os.path.join(root_dir, "docs", "source", "description.txt")
about_workspace_path = os.path.join(root_dir, "docs", "source", "about_workspace.txt")
about_mws_path = os.path.join(root_dir, "docs", "source", "about_mws.txt")
about_plugin_path = os.path.join(root_dir, "docs", "source", "about_plugin.txt")
about_context_path = os.path.join(root_dir, "docs", "source", "about_context.txt")
about_run_path = os.path.join(root_dir, "docs", "source", "about_run.txt")

# Read description files
description, err_msg = open_file(description_path, return_err=True)
about_workspace, err_msg = open_file(about_workspace_path, return_err=True)
about_mws, err_msg = open_file(about_mws_path, return_err=True)
about_plugin, err_msg = open_file(about_plugin_path, return_err=True)
about_context, err_msg = open_file(about_context_path, return_err=True)
about_run, err_msg = open_file(about_run_path, return_err=True)

# If a files can't be read, set the content to an error message
if not description:
    description = f"<error>(openAD description not found)</error>"
else:
    description = description.strip()
if not about_workspace:
    about_workspace = "<error>(Workspace description not found)</error>"
else:
    about_workspace = about_workspace.strip()
if not about_mws:
    about_mws = "<error>(MWS description not found)</error>"
else:
    about_mws = about_mws.strip()
if not about_plugin:
    about_plugin = "<error>(Plugin description not found)</error>"
else:
    about_plugin = about_plugin.strip()
if not about_context:
    about_context = "<error>(Context description not found)</error>"
else:
    about_context = about_context.strip()
if not about_run:
    about_run = "<error>(Run description not found)</error>"
else:
    about_run = about_run.strip()


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
