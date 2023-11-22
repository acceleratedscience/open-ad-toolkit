""" Perform Predict Reaction on a Reaction String"""
from time import sleep
import importlib.util as ilu
from openad.helpers.output import output_table
from openad.helpers.output import output_text
from openad.helpers.output import output_error

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem


def get_include_lib(cmd_pointer):
    """Load RXN include Libraries"""
    folder = cmd_pointer.toolkit_dir + "/RXN" + "/rxn_include.py"
    file = "rxn_include"
    spec = ilu.spec_from_file_location(file, folder)
    rxn = ilu.module_from_spec(spec)
    spec.loader.exec_module(rxn)
    rxn_helper = rxn.rxn_helper()
    return rxn_helper


def get_reaction_from_smiles(
    reaction_smiles: str,
) -> Chem.rdChemReactions.ChemicalReaction:  # pylint: disable=no-member
    """get reaction image from smiles"""
    return AllChem.ReactionFromSmarts(reaction_smiles, useSmiles=True)  # pylint: disable=no-member


def predict_reaction(inputs: dict, cmd_pointer):
    """predict a reaction from a reaction string
    inputs: pyparsing parser object
    cmd_pointer: runtime class"""
    rxn_helper = get_include_lib(cmd_pointer)
    rxn_helper.sync_up_workspace_name(cmd_pointer)
    rxn_helper.get_current_project(cmd_pointer)

    rxn4chemistry_wrapper = cmd_pointer.login_settings["client"][cmd_pointer.login_settings["toolkits"].index("RXN")]
    # Prepare the data query
    rxn_helper = get_include_lib(cmd_pointer)
    molecule = inputs["molecule"]
    val = "val"
    predict_id = None
    ai_model = "2020-08-10"
    if "prediction_id" in inputs:
        predict_id = inputs["prediction_id"][val]
    if "ai_model" in inputs:
        ai_model = inputs["ai_model"][val]

    error_list = []
    if "use_saved" in inputs:
        use_saved = True
    else:
        use_saved = False

    for i in molecule.split("."):
        if not rxn_helper.valid_smiles(str(i)):
            error_list.append(i)

    if len(error_list) > 0:
        df = pd.DataFrame(error_list, columns=["smiles"])
        output_error(" The following invalid Smiles were supplied:", cmd_pointer=cmd_pointer, return_val=False)
        output_table(df, cmd_pointer=cmd_pointer)
        return False
    ##################################################################################################################
    # check for cached reaction
    entry_2 = []

    for i in molecule.split("."):
        entry_2.append(Chem.MolToSmiles(Chem.MolFromSmiles(i), canonical=True))  # pylint: disable=no-member

    result = rxn_helper.retrieve_cache(cmd_pointer, entry_2, "predict_Model-" + ai_model)
    sources = ""
    if result is not False and use_saved is True:
        predict_reaction_results = result
        smiles = predict_reaction_results["response"]["payload"]["attempts"][0]["smiles"]

        confidence = predict_reaction_results["response"]["payload"]["attempts"][0]["confidence"]
        x_y = smiles.split(">>")[1]
        output_text("", cmd_pointer=cmd_pointer, return_val=False)
        output_text("<h2>Saved Result</h2> ", cmd_pointer=cmd_pointer, return_val=False)
        for x in result:
            if len(sources) > 0:
                sources = sources + " + " + x
            else:
                sources = x
            output_text("<success>Smiles:</success>    " + smiles, cmd_pointer=cmd_pointer, return_val=False)
            output_text(
                "<success>Reaction:</success> " + sources + "    ---->    " + x_y,
                cmd_pointer=cmd_pointer,
                return_val=False,
            )
            output_text("<success>Confidence:</success> " + str(confidence), cmd_pointer=cmd_pointer, return_val=False)
        if cmd_pointer.notebook_mode is True:
            return get_reaction_from_smiles(smiles)
        else:
            output_text("", cmd_pointer=cmd_pointer, return_val=False)
            return True

    ##################################################################################################################
    # Since Not Cached lets generate
    try:
        if predict_id is None and ai_model is None:
            predict_reaction_response = rxn4chemistry_wrapper.predict_reaction(molecule)
        elif predict_id is None:
            predict_reaction_response = rxn4chemistry_wrapper.predict_reaction(molecule, ai_model=ai_model)
        else:
            predict_reaction_response = rxn4chemistry_wrapper.predict_reaction(molecule, prediction_id=predict_id)
    except Exception as e:
        raise Exception(
            "\n" + "Prediction call raised an error: " + str(e) + "\n"
        ) from e  # pylint: disable=broad-exception-raised
    status = False
    retries = 0
    while status is False:
        try:
            if predict_id is None and ai_model is None:
                predict_reaction_response = rxn4chemistry_wrapper.predict_reaction(molecule)
            elif predict_id is None:
                predict_reaction_response = rxn4chemistry_wrapper.predict_reaction(molecule, ai_model=ai_model)
            else:
                predict_reaction_response = rxn4chemistry_wrapper.predict_reaction(molecule, prediction_id=predict_id)
            status = True

        except Exception as e:  # pylint: disable=broad-exception-caught
            retries = retries + 1
            sleep(2)
            if retries > 4:
                raise Exception("Server unresponsive" + str(e)) from e  # pylint: disable=broad-exception-raised
    status = False
    sleep(2)
    predict_reaction_results = {}
    retries = 0
    while "response" not in predict_reaction_results:
        try:
            predict_reaction_results = rxn4chemistry_wrapper.get_predict_reaction_results(
                predict_reaction_response["prediction_id"]
            )
            if "response" not in predict_reaction_results:
                sleep(2)
        except Exception as e:  # pylint: disable=broad-exception-caught
            retries = retries + 1
            if retries > 4:
                raise Exception("Server unresponsive" + str(e)) from e  # pylint: disable=broad-exception-raised

    smiles = predict_reaction_results["response"]["payload"]["attempts"][0]["smiles"]

    confidence = predict_reaction_results["response"]["payload"]["attempts"][0]["confidence"]
    rxn_helper.save_to_results_cache(
        cmd_pointer, smiles.split(">>")[0].split("."), predict_reaction_results, "predict_Model-" + ai_model
    )

    source = []
    for i in predict_reaction_results["response"]["payload"]["attempts"][0]["smiles"].split(">>")[0].split("."):
        source.append(i)
    x_y = predict_reaction_results["response"]["payload"]["attempts"][0]["smiles"].split(">>")[1]

    output_text("", cmd_pointer=cmd_pointer, return_val=False)
    output_text("<h2>Generated Result</h2> ", cmd_pointer=cmd_pointer, return_val=False)
    output_text("<success>Smiles:</success>    " + smiles, cmd_pointer=cmd_pointer, return_val=False)
    sources = ""
    for x in source:
        if len(sources) > 0:
            sources = sources + " + " + x
        else:
            sources = x

    output_text(
        "<success>Reaction:</success> " + sources + "    ---->    " + x_y, cmd_pointer=cmd_pointer, return_val=False
    )
    output_text("<success>Confidence:</success> " + str(confidence), cmd_pointer=cmd_pointer, return_val=False)
    if cmd_pointer.notebook_mode is True:
        return get_reaction_from_smiles(predict_reaction_results["response"]["payload"]["attempts"][0]["smiles"])
    else:
        output_text("", cmd_pointer=cmd_pointer, return_val=False)
        return True


def confirm_prompt(question: str) -> bool:
    """confirmation Prompt function"""
    import readline

    reply = None
    while reply not in ("y", "n"):
        reply = input(f"{question} (y/n): ").casefold()
        readline.remove_history_item(readline.get_current_history_length() - 1)
    return reply == "y"
