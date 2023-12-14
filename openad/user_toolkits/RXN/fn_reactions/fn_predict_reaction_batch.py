""" Performs Reaction Prediction on a list of Reactions
"""

from openad.helpers.output import output_table
from openad.helpers.output import output_text
from openad.helpers.output import output_error
from openad.helpers.output import output_warning
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from time import sleep


def get_reaction_from_smiles(reaction_smiles: str) -> Chem.rdChemReactions.ChemicalReaction:
    """gets reactions from a smiles reaction string"""
    return AllChem.ReactionFromSmarts(reaction_smiles, useSmiles=True)  # pylint: disable=no-member


def get_include_lib(cmd_pointer):
    """loads the RXN include module"""
    import importlib.util as ilu

    folder = cmd_pointer.toolkit_dir + "/RXN" + "/rxn_include.py"
    file = "rxn_include"
    spec = ilu.spec_from_file_location(file, folder)
    rxn = ilu.module_from_spec(spec)
    spec.loader.exec_module(rxn)
    rxn_helper = rxn.rxn_helper()
    return rxn_helper


def predict_reaction_batch(inputs: dict, cmd_pointer):
    """Predicts Reactions in Batch from a list of Reaction Smiles Strings"""
    if cmd_pointer.notebook_mode is True:
        from halo import HaloNotebook as Halo  # pylint: disable=import-outside-toplevel
    else:
        from halo import Halo  # pylint: disable=import-outside-toplevel

    rxn_helper = get_include_lib(cmd_pointer)
    rxn_helper.sync_up_workspace_name(cmd_pointer)
    rxn_helper.get_current_project(cmd_pointer)
    if cmd_pointer.notebook_mode is True:
        from IPython.display import display

    class Spinner(Halo):
        """custom spinner"""

        def __init__(self):
            # Alternative spinners:
            # simpleDotsScrolling, interval=100
            super().__init__(spinner="dots", color="white")

    ###################################################################################################
    # getting our input source for the reactions

    if isinstance(inputs["from_source"], dict) and inputs["from_source"]["from_list"] != None:
        try:
            from_list = inputs["from_source"]["from_list"]
        except:  # pylint: disable=bare-except
            output_error(
                "unexpected pyparsing error. Please screenshot and report circumstance to OpenAD team",
                cmd_pointer=cmd_pointer,
                return_val=False,
            )
            output_error("Restart Notebook Kernel or application to proceed", cmd_pointer=cmd_pointer, return_val=False)
            return False

    elif "from_list" in inputs["from_source"][0]:
        try:
            from_list = inputs["from_source"][0]["from_list"]
        except:  # pylint: disable=bare-except
            output_error(
                "unexpected pyparsing error. Please screenshot and report circumstance to OpenAD team",
                cmd_pointer=cmd_pointer,
                return_val=False,
            )
            output_error("Restart Notebook Kernel or application to proceed", cmd_pointer=cmd_pointer, return_val=False)
            return False
    elif "from_dataframe" in inputs:
        try:
            react_frame = cmd_pointer.api_variables[inputs["from_dataframe"]]
            from_list = rxn_helper.get_column_as_list_from_dataframe(react_frame, "reactions")
            if from_list == []:
                raise Exception(
                    "No Provided reactions, data frame should have column 'reactions' "
                )  # pylint: disable=broad-exception-raised
        except Exception:  # pylint: disable=broad-exception-caught
            output_error(
                "Could not load valid list from dataframe column 'reactions' ",
                cmd_pointer=cmd_pointer,
                return_val=False,
            )
            return True
    elif "from_file" in inputs:
        from_file = inputs["from_file"]
        try:
            react_frame = rxn_helper.get_dataframe_from_file(cmd_pointer, from_file)

            from_list = rxn_helper.get_column_as_list_from_dataframe(react_frame, "reactions")
            if from_list == []:
                raise Exception(
                    "No Provided reactions, file should have column 'reactions' "
                )  # pylint: disable=broad-exception-raised
        except Exception:  # pylint: disable=broad-exception-caught
            output_error(
                "Could not load valid list from file column 'reactions' ", cmd_pointer=cmd_pointer, return_val=False
            )
            return True
    newspin = Spinner()

    ### setting up default values... note to put into json metdata file in future

    val = "val"
    if "ai_model" in inputs:
        ai_model = inputs["ai_model"][val]
    else:
        ai_model = "2020-08-10"
    new_from_list = []
    cached_results = []
    if "use_saved" in inputs:
        use_saved = True
    else:
        use_saved = False

    for entry in from_list:
        error_list = []
        for i in entry.split("."):
            if not rxn_helper.valid_smiles(str(i)):
                error_list.append(i)

        if len(error_list) > 0:
            df = pd.DataFrame(error_list, columns=["smiles"])
            output_error(" The following invalid were Smiles Supplied:", cmd_pointer=cmd_pointer, return_val=False)
            output_table(df, cmd_pointer=cmd_pointer)
            output_warning(" This reaction will be skipped  " + entry + " ", cmd_pointer=cmd_pointer, return_val=False)
            continue

        entry_2 = []
        for i in entry.split("."):
            entry_2.append(Chem.MolToSmiles(Chem.MolFromSmiles(i), canonical=True))  # pylint: disable=no-member
        result = rxn_helper.retrieve_cache(cmd_pointer, entry_2, "predict_batch_Model-" + ai_model)

        if result != False and use_saved == True:
            cached_results.append(result)
        else:
            new_from_list.append(entry)

    for reaction_prediction in cached_results:
        source = []
        for i in reaction_prediction["smiles"].split(">>")[0].split("."):
            source.append(i)
        x_y = reaction_prediction["smiles"].split(">>")[1]
        output_text("\n<h2>Saved Result</h2> ", cmd_pointer=cmd_pointer, return_val=False)
        output_text(
            f'<success>Smiles:</success> {reaction_prediction["smiles"]}', cmd_pointer=cmd_pointer, return_val=False
        )
        sources = ""
        for x in source:
            if len(sources) > 0:
                sources = sources + " + " + x
            else:
                sources = x

        output_text(
            f'<success>Confidence:</success> {reaction_prediction["confidence"]}',
            cmd_pointer=cmd_pointer,
            return_val=False,
        )

        if cmd_pointer.notebook_mode is True:
            from IPython.display import display

            display(get_reaction_from_smiles(reaction_prediction["smiles"]))

    if len(new_from_list) > 0:
        newspin.start("Starting Prediction")
        from_list = new_from_list
        rxn4chemistry_wrapper = cmd_pointer.login_settings["client"][
            cmd_pointer.login_settings["toolkits"].index("RXN")
        ]
        retries = 0
        status = False
        while status == False:
            try:
                newspin.text = "Processing Prediction"

                predict_reaction_batch_response = rxn4chemistry_wrapper.predict_reaction_batch(from_list)
                sleep(2)
                status = True
            except Exception as e:  # pylint: disable=broad-exception-caught
                retries = retries + 1
                if retries > 4:
                    newspin.fail("Unable to Process")
                    newspin.stop()
                    raise Exception("Server unresponsive" + str(e)) from e  # pylint: disable=broad-exception-raised

        retries = 0

        reaction_predictions = {}
        while "predictions" not in reaction_predictions:
            try:
                newspin.text = "Processing Prediction"

                reaction_predictions = rxn4chemistry_wrapper.get_predict_reaction_batch_results(
                    predict_reaction_batch_response["task_id"]
                )
                if "predictions" not in reaction_predictions:
                    sleep(3)
            except Exception as e:  # pylint: disable=broad-exception-caught
                retries = retries + 1
                if retries > 10:
                    newspin.fail("Unable to Process")
                    newspin.stop()
                    raise BaseException("Server unresponsive" + str(e)) from e  # pylint: disable=broad-exception-raised
        newspin.succeed("Finished Processing")
        newspin.stop()
        if cmd_pointer.notebook_mode is True:
            from IPython.display import display  # pylint: disable=import-outside-toplevel
        for reaction_prediction in reaction_predictions["predictions"]:
            rxn_helper.save_to_results_cache(
                cmd_pointer,
                reaction_prediction["smiles"].split(">>")[0].split("."),
                reaction_prediction,
                "predict_batch_Model-" + ai_model,
            )
            source = []
            for i in reaction_prediction["smiles"].split(">>")[0].split("."):
                source.append(i)
            x_y = reaction_prediction["smiles"].split(">>")[1]

            output_text("\n<h2>Generated Result</h2> ", cmd_pointer=cmd_pointer, return_val=False)
            output_text(
                f'<success>Smiles:</success> {reaction_prediction["smiles"]}', cmd_pointer=cmd_pointer, return_val=False
            )
            sources = ""
            for x in source:
                if len(sources) > 0:
                    sources = sources + " + " + x
                else:
                    sources = x
            output_text(
                "<success>Reaction:</success> " + sources + "    ---->    " + x_y,
                cmd_pointer=cmd_pointer,
                return_val=False,
            )
            output_text(
                f'<success>Confidence:</success> {reaction_prediction["confidence"]}',
                cmd_pointer=cmd_pointer,
                return_val=False,
            )

            if cmd_pointer.notebook_mode is True:
                display(get_reaction_from_smiles(reaction_prediction["smiles"]))

    output_text(" ", cmd_pointer=cmd_pointer, return_val=False)
    return True
