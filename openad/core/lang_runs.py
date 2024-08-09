"""Contains Run Related functions """

import os
import glob
import readline
from IPython.display import display

# Global variables
from openad.app.global_var_lib import GLOBAL_SETTINGS

# Helpers
from openad.helpers.output import output_text, output_error, output_success, output_table
from openad.helpers.output_msgs import msg


# Create a directory inside the workspace in case it doesn't exit yet.
def _create_workspace_dir_if_nonexistent(cmd_pointer, dir_name):
    if not os.path.isdir(cmd_pointer.workspace_path(cmd_pointer.settings["workspace"].upper()) + "/" + dir_name):
        os.mkdir(cmd_pointer.workspace_path(cmd_pointer.settings["workspace"].upper()) + "/" + dir_name)


# Save a run to the workspace's _runs directory.
def save_run(cmd_pointer, parser):
    """Saves a Run"""
    _create_workspace_dir_if_nonexistent(cmd_pointer, "_runs")
    readline.write_history_file(cmd_pointer.histfile)

    # f =_meta_workspaces+'/'+ cmd_pointer.settings['workspace'].upper()+'/.cmd_history'
    runlist = []

    rows = readline.get_current_history_length() - 1

    while rows > 0 and " ".join(readline.get_history_item(rows).lower().split()) != "create run":
        if " ".join(readline.get_history_item(rows).lower().split()).startswith("save run"):
            return output_error(msg("fail_run_create"))
        runlist.insert(0, readline.get_history_item(rows))
        rows = rows - 1

    f = (
        cmd_pointer.workspace_path(cmd_pointer.settings["workspace"].upper())
        + "/_runs/"
        + parser.as_dict()["run_name"]
        + ".run"
    )
    if rows == 1 and " ".join(readline.get_history_item(rows).lower().split()) != "create run":
        return output_error(msg("fail_run_create"))
    run_file = open(f, "w", encoding="utf-8")

    for i in runlist:
        run_file.write(i + "\n")
    run_file.close()

    return output_success(msg("success_run_save"))


def remove_run(cmd_pointer, parser):
    run_name = parser.as_dict()["run_name"]
    run_file_path = (
        cmd_pointer.workspace_path(cmd_pointer.settings["workspace"].upper()) + "/_runs/" + run_name + ".run"
    )
    if os.path.exists(run_file_path):
        os.remove(run_file_path)
        output_text(f"run {run_name} removed", pad=1)
    else:
        output_text(f"run {run_name} does not exist", pad=1)


# executes a run file.
def exec_run(cmd_pointer, parser):
    """executes a run"""
    _create_workspace_dir_if_nonexistent(cmd_pointer, "_runs")
    f = cmd_pointer.workspace_path(cmd_pointer.settings["workspace"].upper()) + "/.cmd_history"
    f = (
        cmd_pointer.workspace_path(cmd_pointer.settings["workspace"].upper())
        + "/_runs/"
        + parser.as_dict()["run_name"]
        + ".run"
    )
    run_file = open(f, "r", encoding="utf-8")
    run_line = run_file.readline()

    while run_line:
        run_line = run_line.strip()
        if run_line != "":  # Ignore blank lines.
            try:
                # This is to cover a silly edge case when you'd have `?` or `??` in a run.
                # Not  relevant in the real world but used for testing.
                if run_line == "?":
                    output = cmd_pointer.do_help("")
                elif run_line == "??":
                    output = cmd_pointer.default("?")
                else:
                    output = cmd_pointer.default(run_line)

                # When you run a run inside another run, there is no result returned.
                # This prevents "None" from being displayed in Jupyter.
                if GLOBAL_SETTINGS["display"] == "notebook" and output is not None:
                    display(output)
            except Exception as err:
                return output_error(msg("fail_run_execute", run_line, err))

        run_line = run_file.readline()
    run_file.close()


# Display the contents of a Run.
def display_run(cmd_pointer, parser):
    """displays the commands in a run"""
    run_name = parser.as_dict()["run_name"]

    # Create _runs directory if it does not exist yet.
    _create_workspace_dir_if_nonexistent(cmd_pointer, "_runs")

    # import readline
    # readline.write_history_file(cmd_pointer.histfile)  # @Phil, I put this back but it was commented out

    # Read the run file.
    commands = []
    run_file_path = (
        cmd_pointer.workspace_path(cmd_pointer.settings["workspace"].upper()) + "/_runs/" + run_name + ".run"
    )
    try:
        run_file = open(run_file_path, "r", encoding="utf-8")
    except FileNotFoundError:
        return output_error(msg("fail_run_display", run_name))
    line = run_file.readline()
    while line:
        if GLOBAL_SETTINGS["display"] == "notebook":
            commands.append(list([line.replace("\n", "")]))
        else:
            commands.append(line.replace("\n", ""))

        line = run_file.readline()
    run_file.close()
    if GLOBAL_SETTINGS["display"] == "terminal":
        commands = list(
            map(
                lambda t: [output_text(f'<soft>{t[0]}</soft>  { "<cmd>"+str(t[1])+"</cmd>" }', return_val=True)],
                enumerate(commands),
            )
        )

    table_headers = [output_text(f"<soft>Run:</soft> {run_name}", return_val=True)]

    # Display/return table.
    return output_table(commands, is_data=False, headers=table_headers)


# Lists all Runs in the current workspace.
def list_runs(cmd_pointer, parser):
    """list all runs"""
    runs = []

    _create_workspace_dir_if_nonexistent(cmd_pointer, "_runs")

    # Assemble table data.
    for i in glob.glob(
        cmd_pointer.workspace_path(cmd_pointer.settings["workspace"].upper()) + "/_runs/*.run", recursive=True
    ):
        runs.append(list([str(os.path.basename(i).split(".")[0])]))

    # No runs saved yet.
    if len(runs) == 0:
        return output_text(msg("no_runs_saved_yet"), pad=1)

    # Display/return table.
    return output_table(runs, is_data=False, headers=["Stored Runs"])
