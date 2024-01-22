# Example command:
# predict reaction topn in batch from list ['BrBr.c1ccc2cc3ccccc3cc2c1CCO' , 'BrBr.c1ccc2cc3ccccc3cc2c1']

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from time import sleep
import importlib.util as include
from openad.plugins.style_parser import strip_tags
from openad.app.global_var_lib import GLOBAL_SETTINGS
from openad.helpers.spinner import spinner
from openad.helpers.output import output_text, output_warning, output_error, output_table
from openad.helpers.output_msgs import msg
from openad.helpers.general import print_separator
from openad.helpers.general import load_tk_module


def get_reaction_from_smiles(
    reaction_smiles: str,
) -> Chem.rdChemReactions.ChemicalReaction:  # pylint: disable=no-member
    """get a reaction from a smiles defined reaction string"""
    return AllChem.ReactionFromSmarts(reaction_smiles, useSmiles=True)  # pylint: disable=no-member


def predict_reaction_batch_topn(inputs: dict, cmd_pointer):
    """
    Predict top (n) reactions in Batch from a given list of reactions.
    """

    # Load module from toolkit folder
    rxn_helper = load_tk_module(cmd_pointer, "RXN", "rxn_include", "rxn_helper")()

    top_n = 5
    rxn_helper.sync_up_workspace_name(cmd_pointer)
    rxn_helper.get_current_project(cmd_pointer)

    if GLOBAL_SETTINGS["display"] == "notebook":
        from halo import HaloNotebook as Halo  # pylint: disable=import-outside-toplevel
    else:
        from halo import Halo  # pylint: disable=import-outside-toplevel

    val = "val"
    if "topn" in inputs:
        if int(inputs["topn"][val]) > 5:
            top_n = int(inputs["topn"][val])

    if "ai_model" in inputs:
        ai_model = inputs["ai_model"][val]
    else:
        ai_model = "2020-08-10"

    ###################################################################################################
    # getting our input source for the reactions

    if isinstance(inputs["from_source"], dict) and inputs["from_source"]["from_list"] is not None:
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
                    "No Provided reactions, data frame should have column 'reactions'"
                )  # pylint: disable=broad-exception-raised
        except Exception:  # pylint: disable=broad-exception-caught
            output_error(
                "Could not load valid list from dataframe column 'reactions' ",
                return_val=False,
            )
            return False
    elif "from_file" in inputs:
        from_file = inputs["from_file"]
        try:
            react_frame = rxn_helper.get_dataframe_from_file(cmd_pointer, from_file)
            from_list = rxn_helper.get_column_as_list_from_dataframe(react_frame, "reactions")
            if from_list == []:
                raise Exception  # pylint: disable=broad-exception-raised
        except Exception:  # pylint: disable=broad-exception-caught
            output_error("Could not load valid list from file column 'reactions' ", return_val=False)
            return False

    if "use_saved" in inputs:
        use_saved = True
    else:
        use_saved = False

    new_from_list = []
    cached_results = []
    new_cannonical_list = []
    cached_cannonical_list = []

    for entry in from_list:
        error_list = []
        for i in entry.split("."):
            if not rxn_helper.valid_smiles(str(i)):
                error_list.append(i)

        if len(error_list) > 0:
            df = pd.DataFrame(error_list, columns=["smiles"])
            output_error(":The following invalid were Smiles Supplied:", return_val=False)
            output_table(df)
            output_warning("info: This rection will be skipped: " + entry + " ", return_val=False)
            continue

        entry_2 = []
        for i in entry.split("."):
            entry_2.append(Chem.MolToSmiles(Chem.MolFromSmiles(i), canonical=True))  # pylint: disable=no-member

        result = rxn_helper.retrieve_cache(
            cmd_pointer, entry_2, "predict_batch_topn" + str(top_n) + "_model" + ai_model
        )
        if result is not False and use_saved is True:
            cached_results.append(result)
            cached_cannonical_list.append(entry_2)
        else:
            new_from_list.append(entry.split("."))
            new_cannonical_list.append(entry_2)
    reaction_no = 0
    output_text("\n", return_val=False)
    # address already saved results
    for reaction_predictions in cached_results:
        output_text(
            f"<green>Saved Reaction results for:</green> {cached_cannonical_list[reaction_no]}",
            return_val=False,
        )

        reaction_no = reaction_no + 1
        for j, prediction in enumerate(reaction_predictions["results"], 1):
            product_smiles = ".".join(prediction["smiles"])
            confidence = prediction["confidence"]
            output_text(
                f"<green>         Product(s) </green>{j} {product_smiles}, With confidence {confidence}",
                return_val=False,
            )
        output_text("\n", return_val=False)
    # if requirement for ask RXN to generate Results
    if len(new_from_list) > 0:
        val = "val"
        from_list = new_from_list
        retries = 0
        status = False
        rxn4chemistry_wrapper = cmd_pointer.login_settings["client"][
            cmd_pointer.login_settings["toolkits"].index("RXN")
        ]

        spinner.start("Processing")

        while status is False:
            try:
                # spinner.text = "Processing Prediction"
                sleep(2)
                predict_rection_batch_response = rxn4chemistry_wrapper.predict_reaction_batch_topn(
                    precursors_lists=new_from_list,
                    topn=top_n,
                    ai_model=ai_model,
                )
                status = True
            except Exception as err:  # pylint: disable=broad-exception-caught
                retries = retries + 1
                if retries > 4:
                    spinner.fail("Unable to process input")
                    spinner.stop()
                    output_error(["Server unresponsive", err], return_val=False)
                    return

        x = {}
        retries = 0
        while "predictions" not in x:
            try:
                # spinner.text = "Processing Prediction"
                x = rxn4chemistry_wrapper.get_predict_reaction_batch_topn_results(
                    predict_rection_batch_response["task_id"]
                )
                if "predictions" not in x:
                    sleep(3)
            except Exception as e:  # pylint: disable=broad-exception-caught
                retries = retries + 1
                if retries > 10:
                    spinner.fail("Unable to Process")
                    spinner.stop()
                    raise Exception("Server unresponsive " + str(e)) from e  # pylint: disable=broad-exception-raised

        reaction_no = 0
        spinner.stop()
        output_final = []
        alphabet = [
            "A",
            "B",
            "C",
            "D",
            "E",
            "F",
            "G",
            "H",
            "I",
            "J",
            "K",
            "L",
            "M",
            "N",
            "O",
            "P",
            "Q",
            "R",
            "S",
            "T",
            "U",
            "V",
            "W",
            "X",
            "Y",
            "Z",
            "AA",
            "BB",
            "CC",
            "DD",
            "EE",
            "FF",
            "GG",
            "HH",
            "II",
            "JJ",
            "KK",
            "LL",
            "MM",
            "NN",
            "OO",
            "PP",
            "QQ",
            "RR",
            "SS",
            "TT",
            "UU",
            "VV",
            "WW",
            "XX",
            "YY",
            "ZZ",
        ]

        # NOT USED BUT MAYBE LATER - see below
        # table = []
        # table_headers=['R', '#', 'Input A', 'Input B', 'Prediction', 'Confidence']

        for i, reaction_predictions in enumerate(x["predictions"], 1):
            input = new_cannonical_list[reaction_no]
            output = []

            # NOT USED BUT MAYBE LATER - see below
            # row_shared = [alphabet[i - 1], None, input, input[1]] # Would have to be updated to support > 2

            output.append(f"Reaction predictions for:")
            for i, mol in enumerate(input):
                output.append(f"{alphabet[i]}: {mol}")
            output.append("---")
            rxn_helper.save_to_results_cache(
                cmd_pointer,
                input,
                reaction_predictions,
                "predict_batch_topn" + str(top_n) + "_model" + ai_model,
            )
            reaction_no = reaction_no + 1

            for j, prediction in enumerate(reaction_predictions["results"], 1):
                product_smiles = ".".join(prediction["smiles"])
                confidence = prediction["confidence"]
                confidence = round(confidence * 100)

                if confidence < 1:
                    confidence_str = f"<1"
                else:
                    confidence_str = f"{confidence:2d}"

                output_line = f"{j} - confidence {confidence_str}% - {product_smiles}"

                # Colorize list results according to confidence.
                if confidence > 66:
                    output_line = f"<green>{output_line}</green>"
                elif confidence > 33:
                    output_line = f"<yellow>{output_line}</yellow>"
                elif confidence < 10 > 1:
                    output_line = f"<soft>{output_line}</soft>"
                elif confidence < 1:
                    output_line = f"<soft>{output_line}</soft>"

                output.append(output_line)

                # NOT USED BUT MAYBE LATER - see below
                # # Create table row
                # row = row_shared.copy()
                # row[1] = j
                # row.append(product_smiles)
                # row.append(confidence)
                #
                # # Colorize table results according to confidence.
                # if confidence > 66:
                #     for idx in range(len(row)):
                #         row[idx] = f'<green>{str(row[idx])}</green>'
                # elif confidence > 33:
                #     for idx in range(len(row)):
                #         row[idx] = f'<yellow>{str(row[idx])}</yellow>'
                # elif confidence < 10 > 1:
                #     for idx in range(len(row)):
                #         row[idx] = f'<soft>{str(row[idx])}</soft>'
                # elif confidence < 1:
                #     for idx in range(len(row)):
                #         row[idx] = f'<soft>< 1</soft>'
                # table.append(row)

            # Add separator line.
            list_width = 0
            for line in output:
                ll = len(strip_tags(line))
                if ll > list_width:
                    list_width = ll
            sep = print_separator("yellow", list_width, return_val=True)
            output[len(input) + 1] = sep

            # Add reaction to main output.
            br = ["", ""] if i > 1 else []
            output_final = output_final + br + output

        output_text("\n".join(output_final), return_val=False, pad_btm=1)

        # Same information as a table. Decided it's not useful, but could possibly
        # be re-enabled via a follow-up command. Something like `result as table`
        # output_table(table, headers=table_headers, return_val=False)
    return True
