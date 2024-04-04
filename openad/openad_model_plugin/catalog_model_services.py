import pyparsing as py
import os
import glob
from openad.helpers.output import output_text, output_table, output_warning, output_error, output_success
SERVICE_DEFINTION_PATH = '/definition/services/'


def help_dict_create(
    name: str,  # Name of the comand - used for ...?
    command: str,  # Command structure, used for help, docs, training
    description: str,  # Description of the command, used for help, docs, training
    note: str = None,  # Additional note to the command, only used in help (eg. To learn more about runs, run `run ?`)
    url: str = None,  # Currently not used - URL to the documentation of the command?
    category: str = "Uncategorized",  # Category used to organize the commands in help & docs
    parent: str = None,  # Parent command, only relevant for follow-up commands like `result open`
):
    """Create a help dictionary"""
    return {
        "category": category,
        "name": name,
        "command": command,
        "description": description,
        "note": note,
        "url": url,
        "parent": parent,
    }


def list_cataloged_services(cmd_pointer, parser):
    pass


def catalog_service(cmd_pointer, parser):
    instruction = parser.to_dict()
    if "service_name"
    if "path" in instruction:
        if os.path.exists(instruction["path"]+SERVICE_DEFINTION_PATH):

        else:
            return False
    else:
        return False


    pass


def service_up(cmd_pointer, parser):
    pass


def service_down(cmd_pointer, parser):
    pass


def remove_cataloged_service(cmd_pointer, parser):
    pass


def service_catalog_grammar(statements: list, help: list, service_list: list):

    add = py.CaselessKeyword("add")
    model = py.CaselessKeyword("model")
    service = py.CaselessKeyword("service")
    fr_om = py.CaselessKeyword("from")
    path = py.CaselessKeyword("path")
    quoted_string = py.QuotedString("'", escQuote="\\")
    a_s = py.CaselessKeyword("as")

    statements.append(py.Forward(add + model + service + fr_om + path+quoted_string("path")+a_s+quoted_string('service_name'))("add_model_path"))
    help.append(
        help_dict_create(
            name="Add Model from Path",
            category="General",
            command="add model from path '<path to model directory>' as '<service_name>'",
            description="add a model definition to the catalog.",
        )
    )
