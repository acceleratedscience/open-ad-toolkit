"""Welcome messages for each module."""


adccl_intro: str = """

<h1>Introduction to the ADCCL</h1>

ADCCL (Accelerated Discovery Common Client) is a CLI providing centralized access to a number of open-source scientific toolkits developed by IBM Research. These toolkits let you combine the power of artificial intelligence, computational chemistry and cloud computing to supercharge your discovery process and drastically accelerate your development timelines.

<link>http://accelerate.science</link>


<h2>Workspace</h2>
    A workspace is a local sandbox environment where you can experiment at will. You can set up as many workspaces as you need, and each workspace comes with its own settings and history.

    To learn more:
        <cmd>workspace ?</cmd>

<h2>Toolkits</h2>
    ADCCL supports 4 different toolkits as of now.
    • DS4SD - Deep Search
    • GT4SD - Generative Toolkit
    • ST4SD - Simulation Toolkit
    • RXN - Digital Chemistry

    You can access any toolkit's functionality by prepending its name:
        <cmd>DS4SD <do_xyz></cmd>  <green># To be implemented</green>

    To learn more about any tool:
        <cmd>DS4SD</cmd>
        <cmd>GT4SD</cmd>
        <cmd>ST4SD</cmd>
        <cmd>RXN</cmd>

<h2>Contexts</h2>
    By default, you'll be working in the 'root' context. If you are working primarily with one tool, you can switch your context to the tool in question:
        <cmd>set context DS4SD</cmd>
        <cmd><do-xyz></cmd>

    To learn more:
        <cmd>context ?</cmd>

<h2>Getting Started</h2>
    Start with one of our example processes to get started in no time.
    • <link>This is an example link placeholder</link>
"""
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
# adccl_intro = style(adccl_intro, pad=2, tabs=1, edge=True)
