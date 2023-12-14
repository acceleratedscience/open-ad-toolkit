"""Welcome messages for each module."""

info_workspaces = "A workspace is an isolated directory where all of your files, settings and runs are stored. This allows you to work on unrelated projects without them contaminating one another. You can create as many workspaces as you need, and each workspace comes with its own command history."
info_toolkits = """Toolkits are individual pieces of software that provide highly specialized functionality. You can think of OpenAD as the operating system, and toolkits as the applications. Toolkits can be installed and removed independently of one another, and you can switch between them at any time.

As of now, OpenAD supports 4 different toolkits.
• DS4SD - Deep Search
• GT4SD - Generative Toolkit
• ST4SD - Simulation Toolkit
• RXN - Digital Chemistry

To learn more about what any of these toolkits does, run its name: <cmd>ds4sd</cmd> / <cmd>gt4sd</cmd> / <cmd>st4sd</cmd> / <cmd>rxn</cmd>"""
# # Add this back once this functionality is built:
# You can access any toolkit's functionality by prepending its name:
# <cmd>DS4SD <do_xyz></cmd>  <green># To be implemented</green>

info_runs = "A run is a collection of commands that you can save and execute later. Runs are stored in your workspace, and can be shared with other users."

info_context = "After you've installed a toolkit, you need to set your context to that toolkit in order to access its functionality."

openad_intro: str = f"""<h1>Introduction to OpenAD</h1>

OpenAD is a CLI providing centralized access to a number of open-source scientific toolkits developed by IBM Research. These toolkits let you combine the power of artificial intelligence, computational chemistry and hybrid cloud to supercharge your discovery process and drastically accelerate your development timelines.

<link>http://accelerate.science</link>


<h2>Workspaces</h2>
{info_workspaces}

<soft>To see how to work with workspaces, run <cmd>workspace ?</cmd></soft>


<h2>Toolkits</h2>
{info_toolkits}

<soft>To see how to work with toolkits, run <cmd>toolkit ?</cmd></soft>


<h2>Contexts</h2>
{info_context}

For example:
<cmd>add toolkit DS4SD</cmd>
<cmd>set context DS4SD</cmd>
<cmd><do-xyz></cmd>

<h2>Runs</h2>
{info_runs}

<soft>To see how to work with runs, run <cmd>run ?</cmd></soft>"""

# # Add this back once we have some example processes:
# <h2>Getting Started</h2>
# Start with one of our example processes to get started in no time.
# • <link>This is an example link placeholder</link>

# # Project definition was removed because we're not building this yet.
# <h2>Project</h2>
# A project is a cloud environment where you can share models, datasets
# and other files with team members. Every workspace can be associated
# with maximum one project. You can push files from your workspace to a
# project, or pull files from a project into your workspace. Projects have
# version control, so you always know who contributed what, when.

# To learn more:
#     <cmd>project ?</cmd>
#     <cmd>pull ?</cmd>
#     <cmd>push ?</cmd>
# openad_intro = style(openad_intro, pad=2, tabs=1, edge=True)
