""" Interprets a free text paragraph into a set of Recipe instructions """
import importlib.util as ilu
from ad4e_opentoolkit.helpers.output import output_text

list_of_reactions = []


def get_include_lib(cmd_pointer):
    """Import RXN Include Libraries"""
    folder = cmd_pointer.toolkit_dir + "/RXN" + "/rxn_include.py"
    file = "rxn_include"
    spec = ilu.spec_from_file_location(file, folder)
    rxn = ilu.module_from_spec(spec)
    spec.loader.exec_module(rxn)
    rxn_helper = rxn.rxn_helper
    return rxn_helper


def interpret_recipe(inputs: dict, cmd_pointer):
    """Interprets a free text paragraph into a set of Recipe instructions"""
    receipe = inputs["receipe"]
    rxn4chemistry_wrapper = cmd_pointer.login_settings["client"][cmd_pointer.login_settings["toolkits"].index("RXN")]
    # Prepare the data query

    if cmd_pointer.notebook_mode is True:
        from IPython.display import Markdown  # pylint: disable=import-outside-toplevel
    try:
        return_result = []
        actios_from_procedure_results = rxn4chemistry_wrapper.paragraph_to_actions(receipe)
        if cmd_pointer.notebook_mode is True:
            return_result.append("***See the following actions from the Receipe:***\n")
        else:
            return_result.append("See the following actions from the Receipe:\n")

        for index, action in enumerate(actios_from_procedure_results["actions"], 1):
            if cmd_pointer.notebook_mode is True:
                return_result.append(f"{index}. {action}\n")
            else:
                return_result.append(f"{index}. {action}\n")
        if cmd_pointer.notebook_mode is True:
            return Markdown("".join(return_result))
        else:
            return output_text("\n" + "".join(return_result), cmd_pointer=cmd_pointer)
    except Exception as e:
        raise Exception(
            "unable to to turn paragraph to actions:" + str(e)
        ) from e  # pylint: disable=broad-exception-raised
