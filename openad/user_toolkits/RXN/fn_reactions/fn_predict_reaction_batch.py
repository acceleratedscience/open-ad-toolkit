# Example commands:
# predict reaction in batch from list ['BrBr.c1ccc2cc3ccccc3cc2c1CCO' , 'BrBr.c1ccc2cc3ccccc3cc2c1']
# predict reaction in batch from list ['BrBr.c1ccc2cc3ccccc3cc2c1CCO' , 'BrBr.c1ccc2cc3ccccc3cc2c1'] use_saved

"""Performs Reaction Prediction on a list of Reactions"""

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from time import sleep
from openad.app.global_var_lib import GLOBAL_SETTINGS
from openad.helpers.output import output_text, output_error, output_warning, output_table
from openad.helpers.output_msgs import msg
from openad.helpers.general import load_tk_module
from openad.helpers.spinner import Spinner


def get_reaction_from_smiles(reaction_smiles: str) -> Chem.rdChemReactions.ChemicalReaction:
    """gets reactions from a smiles reaction string"""
    return AllChem.ReactionFromSmarts(reaction_smiles, useSmiles=True)  # pylint: disable=no-member


def predict_reaction_batch(inputs: dict, cmd_pointer):
    """Predicts Reactions in Batch from a list of Reaction Smiles Strings"""

    # Load module from toolkit folder
    rxn_helper = load_tk_module(cmd_pointer, "RXN", "rxn_include", "rxn_helper")()

    rxn_helper.sync_up_workspace_name(cmd_pointer)
    rxn_helper.get_current_project(cmd_pointer)
    if GLOBAL_SETTINGS["display"] == "notebook":
        from IPython.display import display

    ###################################################################################################
    # getting our input source for the reactions

    if isinstance(inputs["from_source"], dict) and inputs["from_source"]["from_list"] != None:
        try:
            from_list = inputs["from_source"]["from_list"]
        except Exception:  # pylint: disable=broad-except
            output_error(msg("err_pyparsing"), return_val=False)
            return False

    elif "from_list" in inputs["from_source"][0]:
        try:
            from_list = inputs["from_source"][0]["from_list"]
        except Exception:  # pylint: disable=broad-except
            output_error(msg("err_pyparsing"), return_val=False)
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
            output_error("Could not load valid list from file column 'reactions' ", return_val=False)
            return True
    newspin = Spinner(GLOBAL_SETTINGS["VERBOSE"])

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
            output_error(" The following invalid were Smiles Supplied:", return_val=False)
            output_table(df, is_data=False)
            output_warning(" This reaction will be skipped  " + entry + " ", return_val=False)
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
        output_text("\n<h2>Saved Result</h2> ", return_val=False)
        output_text(f'<green>Smiles:</green> {reaction_prediction["smiles"]}', return_val=False)
        sources = ""
        for x in source:
            if len(sources) > 0:
                sources = sources + " + " + x
            else:
                sources = x

        output_text(
            f'<green>Confidence:</green> {reaction_prediction["confidence"]}',
            return_val=False,
        )

        if GLOBAL_SETTINGS["display"] == "notebook":
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
                if retries == 0:
                    newspin.info("Processing Prediction")

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
        if GLOBAL_SETTINGS["display"] == "notebook":
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

            output_text("\n<h2>Generated Result</h2> ", return_val=False)
            output_text(f'<green>Smiles:</green> {reaction_prediction["smiles"]}', return_val=False)
            sources = ""
            for x in source:
                if len(sources) > 0:
                    sources = sources + " + " + x
                else:
                    sources = x
            output_text(
                "<green>Reaction:</green> " + sources + "    ---->    " + x_y,
                return_val=False,
            )
            output_text(
                f'<green>Confidence:</green> {reaction_prediction["confidence"]}',
                return_val=False,
            )

            if GLOBAL_SETTINGS["display"] == "notebook":
                display(get_reaction_from_smiles(reaction_prediction["smiles"]))

    output_text(" ", return_val=False)
    if not GLOBAL_SETTINGS["VERBOSE"]:
        return reaction_predictions["predictions"]
    else:
        return True
