
import os

# Global variables
from global_var_lib import _meta_dir as _meta_dir
from global_var_lib import _meta_dir_toolkits as _meta_dir_toolkits
from global_var_lib import _meta_registry as _meta_registry
from global_var_lib import _meta_registry_session as _meta_registry_session
from global_var_lib import _meta_login_registry as _meta_login_registry
from global_var_lib import _meta_workspaces as _meta_workspaces
from global_var_lib import _meta_registry_settings as _meta_registry_settings

# Helpers
from helpers.output import msg, output_text, output_error, output_success, output_table
from helpers.general import is_notebook_mode


# Create a directory inside the workspace in case it doesn't exit yet.
def _create_workspace_dir_if_nonexistent(cmd_pointer, dir_name):
    if not os.path.isdir(cmd_pointer.workspace_path(cmd_pointer.settings['workspace'].upper()) + '/' + dir_name):
        os.mkdir(cmd_pointer.workspace_path(cmd_pointer.settings['workspace'].upper()) + '/' + dir_name)


# Save a run to the workspace's _runs directory.
def save_run(cmd_pointer, parser):
    import readline

    _create_workspace_dir_if_nonexistent(cmd_pointer, '_runs')
    readline.write_history_file(cmd_pointer.histfile)

    # f =_meta_workspaces+'/'+ cmd_pointer.settings['workspace'].upper()+'/.cmd_history'
    runlist = []

    rows = readline.get_current_history_length() - 1

    while rows > 0 and ' '.join(readline.get_history_item(rows).lower().split()) != 'create run':
        if ' '.join(readline.get_history_item(rows).lower().split()).startswith('save run'):
            return output_error(msg('fail_run_create'), cmd_pointer)
        runlist.insert(0, readline.get_history_item(rows))
        rows = rows - 1

    f = cmd_pointer.workspace_path(cmd_pointer.settings['workspace'].upper()) + '/_runs/' + parser.as_dict()['run_name'] + '.run'
    if rows == 1 and ' '.join(readline.get_history_item(rows).lower().split()) != 'create run':
        return output_error(msg('fail_run_create'), cmd_pointer)
    run_file = open(f, 'w')
    for i in runlist:
        run_file.write(i + '\n')
    run_file.close()
    return output_success(msg('success_run_save'), cmd_pointer)


# executes a run file.
def execute_run(cmd_pointer, parser):
    _create_workspace_dir_if_nonexistent(cmd_pointer, '_runs')
    f = cmd_pointer.workspace_path(cmd_pointer.settings['workspace'].upper()) + '/.cmd_history'
    f = cmd_pointer.workspace_path(cmd_pointer.settings['workspace'].upper()) + '/_runs/' + parser.as_dict()['run_name'] + '.run'
    run_file = open(f, 'r')
    run_line = run_file.readline()

    while run_line:
        try:
            if cmd_pointer.notebook_mode:
                from IPython.display import display
                display(cmd_pointer.default(run_line))
            else:
                cmd_pointer.default(run_line)
        except BaseException:
            return output_error(msg('fail_run_execute', run_line), cmd_pointer)

        run_line = run_file.readline()
    run_file.close()


# Display the contents of a Run.
def display_run(cmd_pointer, parser):
    table_headers = (f'<soft>Run:</soft> {parser.as_dict()["run_name"]}',)  # NOT Sure why we mess things up with
    notebook_mode = cmd_pointer.notebook_mode
    run_name = parser.as_dict()['run_name']

    # Create _runs directory if it does not exist yet.
    _create_workspace_dir_if_nonexistent(cmd_pointer, '_runs')

    # import readline
    # readline.write_history_file(cmd_pointer.histfile)  # @Phil, I put this back but it was commented out

    # Read the run file.
    commands = []
    run_file_path = cmd_pointer.workspace_path(cmd_pointer.settings['workspace'].upper()) + '/_runs/' + run_name + '.run'
    try:
        run_file = open(run_file_path, 'r')
    except FileNotFoundError:
        return output_error(msg('fail_run_display', run_name, split=True), cmd_pointer)
    line = run_file.readline()
    while line:
        if notebook_mode == True:
            commands.append(list([line.replace('\n', '')]))
        else:
            commands.append(line.replace('\n', ''))

        line = run_file.readline()
    run_file.close()
    if not notebook_mode and not cmd_pointer.api_mode:
        commands = list(
            map(lambda t: [output_text(f'<soft>{t[0]}</soft>  { "<cmd>"+str(t[1])+"</cmd>" }', cmd_pointer, return_val=True)], enumerate(commands)))

    table_headers = [f'<soft>Run:</soft> {run_name}']

    # Display/return table.
    return output_table(commands, cmd_pointer, headers=table_headers)


# Lists all Runs in the current workspace.
def list_runs(cmd_pointer, parser):
    import glob
    runs = []
    notebook_mode = cmd_pointer.notebook_mode
    api_mode = cmd_pointer.api_mode

    _create_workspace_dir_if_nonexistent(cmd_pointer, '_runs')

    # Assemble table data.
    for i in glob.glob(cmd_pointer.workspace_path(cmd_pointer.settings['workspace'].upper()) + '/_runs/*.run', recursive=True):
        runs.append(list([str(os.path.basename(i).split('.')[0])]))

    # No runs saved yet.
    if len(runs) == 0:
        return output_text(msg('no_runs_saved_yet'), cmd_pointer, pad=1)

    # Display/return table.
    return output_table(runs, cmd_pointer, headers=['Stored Runs'])
