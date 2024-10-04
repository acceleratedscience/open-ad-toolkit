# Example command:
# predict retrosynthesis 'BrCCc1cccc2c(Br)c3ccccc3cc12' using (max_steps=3)

"""This library is for implementing the predict retrosynthesis function from RXN"""

from typing import Dict, List
import importlib.util as ilu
from openad.helpers.output import output_text, output_error
from rdkit import Chem
from rdkit.Chem import AllChem
from time import sleep
from openad.app.global_var_lib import GLOBAL_SETTINGS
from openad.smols.smol_cache import create_analysis_record, save_result
from openad.smols.smol_functions import canonicalize, valid_smiles
from openad.helpers.general import load_tk_module
from openad.helpers.spinner import Spinner


def get_reaction_from_smiles(reaction_smiles: str) -> Chem.rdChemReactions.ChemicalReaction:
    """get reaction image"""
    return AllChem.ReactionFromSmarts(reaction_smiles, useSmiles=True)  # pylint: disable=no-member


def collect_reactions_from_retrosynthesis(tree: Dict) -> List[str]:
    """collect all reactions from retrosynthesis tree"""
    reactions = []
    if "children" in tree and len(tree["children"]):
        reactions.append(
            AllChem.ReactionFromSmarts(  # pylint: disable=no-member
                "{}>>{}".format(".".join([node["smiles"] for node in tree["children"]]), tree["smiles"]), useSmiles=True
            )
        )  # pylint: disable=no-member
    for node in tree["children"]:
        reactions.extend(collect_reactions_from_retrosynthesis(node))
    return reactions


def collect_reactions_from_retrosynthesis_text(tree: Dict) -> List[str]:
    """collect all reactions from retrosynthesis tree"""
    reactions = []
    # print(tree)
    if "children" in tree and len(tree["children"]):
        reactions.append(
            "{} --->> {}".format(" + ".join([node["smiles"] for node in tree["children"]]), tree["smiles"])
        )

    for node in tree["children"]:
        reactions.extend(collect_reactions_from_retrosynthesis_text(node))

    return reactions


def predict_retro(inputs: dict, cmd_pointer):
    """Perform RXN Predict Retro Synthesis"""

    # Load module from toolkit folder
    rxn_helper = load_tk_module(cmd_pointer, "RXN", "rxn_include", "rxn_helper")()

    rxn_helper.sync_up_workspace_name(cmd_pointer)
    rxn_helper.get_current_project(cmd_pointer)
    if GLOBAL_SETTINGS["display"] == "notebook":
        from IPython.display import display  # pylint: disable=import-outside-toplevel
        from halo import HaloNotebook as Halo  # pylint: disable=import-outside-toplevel
    else:
        from halo import Halo  # pylint: disable=import-outside-toplevel

    val = "val"
    availability_pricing_threshold = 0
    available_smiles = None
    exclude_smiles = None
    exclude_substructures = None
    exclude_target_molecule = True
    fap = 0.6
    max_steps = 5
    nbeams = 10
    pruning_steps = 2
    ai_model = "2020-07-01"

    product_smiles = inputs["molecule"]
    #######################

    if not valid_smiles(str(product_smiles)):
        output_error(" Invalid Smiles Supplied.", return_val=False)
        return False
    else:
        product_smiles = canonicalize(product_smiles)

    if len(product_smiles.split(".")) > 1:
        output_error(
            " SMILES provides describes a reaction. Use `predict reaction` to see probable result",
            return_val=False,
        )
        return False

    if GLOBAL_SETTINGS["display"] == "notebook":
        import py3Dmol

        style = "stick"
        mol = Chem.MolFromSmiles(product_smiles)  # pylint: disable=no-member
        mol = Chem.AddHs(mol)  # pylint: disable=no-member
        AllChem.EmbedMolecule(mol)  # pylint: disable=no-member
        AllChem.MMFFOptimizeMolecule(mol, maxIters=200)  # pylint: disable=no-member
        mblock = Chem.MolToMolBlock(mol)  # pylint: disable=no-member

        view = py3Dmol.view(width=700, height=500)
        view.addModel(mblock, "mol")
        view.setStyle({style: {}})
        view.zoomTo()
        view.show()
        output_text("<green>Target Molecule:</green> " + product_smiles, return_val=False)
    #######################
    result_parameters = {}
    if "availability_pricing_threshold" in inputs:
        availability_pricing_threshold = int(inputs["availability_pricing_threshold"][val])
        result_parameters["availability_pricing_threshold"] = availability_pricing_threshold
    if "available_smiles" in inputs:
        available_smiles = inputs["available_smiles"][val]
        result_parameters["available_smiles"] = available_smiles
    if "exclude_smiles" in inputs:
        exclude_smiles = inputs["exclude_smiles"][val]
        result_parameters["exclude_smiles"] = exclude_smiles
    if "exclude_substructures" in inputs:
        exclude_substructures = inputs["exclude_substructures"][val]
        result_parameters["exclude_substructures"] = exclude_substructures
    if "exclude_target_molecule" in inputs:
        if inputs["exclude_substructures"][val].upper() == "TRUE":
            exclude_target_molecule = True
            result_parameters["exclude_target_molecule"] = exclude_target_molecule
    if "fap" in inputs:
        fap = float(inputs["fap"][val])
        result_parameters["fap"] = fap
    if "max_steps" in inputs:
        max_steps = int(inputs["max_steps"][val])
        result_parameters["max_steps"] = max_steps
    if "nbeams" in inputs:
        nbeams = int(inputs["nbeams"][val])
        result_parameters["nbeams"] = nbeams
    if "pruning_steps" in inputs:
        pruning_steps = int(inputs["pruning_steps"][val])
        result_parameters["pruning_steps"] = pruning_steps
    if "ai_model" in inputs:
        ai_model = inputs["ai_model"][val]
        result_parameters["ai_model"] = ai_model
    rxn4chemistry_wrapper = cmd_pointer.login_settings["client"][cmd_pointer.login_settings["toolkits"].index("RXN")]

    # Prepare the data query
    print("\n")
    newspin = Spinner(GLOBAL_SETTINGS["VERBOSE"])
    newspin.start("Starting Retrosynthesis")
    try:
        retries = 0
        status = False
        predict_retro_response = None

        while retries < 10 and status is False:
            try:
                if retries == 0:
                    newspin.info("Submitting Retrosynthesis ")
                predict_retro_response = rxn4chemistry_wrapper.predict_automatic_retrosynthesis(
                    product_smiles,
                    availability_pricing_threshold=availability_pricing_threshold,
                    available_smiles=available_smiles,
                    exclude_smiles=exclude_smiles,
                    exclude_substructures=exclude_substructures,
                    exclude_target_molecule=exclude_target_molecule,
                    fap=fap,
                    max_steps=max_steps,
                    nbeams=nbeams,
                    pruning_steps=pruning_steps,
                    ai_model=ai_model,
                )
                status = True
            except Exception as e:  # pylint: disable=broad-exception-caught
                sleep(2)
                retries = retries + 1
                if retries >= 10:
                    newspin.fail("Unable to Process")
                    newspin.start()
                    newspin.stop()
                    raise Exception(
                        "Server unresponsive: Unable to submit for processing after 10 retires" + str(e)
                    ) from e  # pylint: disable=broad-exception-raised

        if predict_retro_response is None:
            raise BaseException(
                "Server unresponsive: Unable to submit for processing after 10 retires"
            )  # pylint: disable=broad-exception-raised

        retries = 0
        if predict_retro_response["response"]["payload"]["errorMessage"] is not None:
            return predict_retro_response["response"]["payload"]["errorMessage"]

        status = "NEW"
        previous_status = status
        while status != "SUCCESS":
            try:
                if status != previous_status:
                    newspin.info("Processing Retrosynthesis :" + status)
                    previous_status = status

                predict_automatic_retrosynthesis_results = (
                    rxn4chemistry_wrapper.get_predict_automatic_retrosynthesis_results(
                        predict_retro_response["prediction_id"]
                    )
                )
                if predict_retro_response["response"]["payload"] is None:
                    output_text(
                        "<h2>No Result:</h2>  Unable to find path for  " + product_smiles,
                        return_val=False,
                    )
                    if GLOBAL_SETTINGS["display"] == "notebook":
                        return
                    else:
                        return False
                status = predict_automatic_retrosynthesis_results["status"]

                sleep(5)

            except Exception as e:  # pylint: disable=broad-exception-caught
                retries = retries + 1
                sleep(15)
                status = "Waiting"
                newspin.info("Processing Retrosynthesis: Waiting")
                if retries > 20:
                    raise Exception(
                        "Server unresponsive: Unable to complete processing for prediction id:'"
                        + predict_retro_response["prediction_id"]
                        + "'after 20 retires"
                        + str(e)
                    ) from e  # pylint: disable=broad-exception-raised
    except Exception as e:
        newspin.fail("Unable to Process")
        newspin.start()
        newspin.stop()
        raise Exception("Unable to complete processing " + str(e)) from e  # pylint: disable=broad-exception-raised
    reactions_text = []
    # print(predict_automatic_retrosynthesis_results)
    try:
        for index, tree in enumerate(predict_automatic_retrosynthesis_results["retrosynthetic_paths"]):
            # print("outer")
            for reaction in collect_reactions_from_retrosynthesis_text(tree):
                reactions_text.append(str(reaction))
                # print("inner")

    except Exception as e:  # pylint: disable=broad-exception-caught
        newspin.fail("Unable to Process")
        newspin.stop()
        raise Exception(
            "The following Error message was received while trying to process results:" + str(e)
        ) from e  # pylint: disable=broad-exception-raised
    num_results = 0
    # print(reactions_text)
    try:
        newspin.succeed("Finished Processing")
        newspin.start()
        newspin.stop()
        results = {}
        i = 0
        for index, tree in enumerate(predict_automatic_retrosynthesis_results["retrosynthetic_paths"]):
            num_results = num_results + 1
            output_text(
                "",
                return_val=False,
            )
            if num_results < 4 or GLOBAL_SETTINGS["VERBOSE"] == False:
                results[str(index)] = {"confidence": tree["confidence"], "reactions": []}

            output_text(
                "<h2> <green>Showing path </green> {} <green> with confidence </green>{}:".format(
                    index, tree["confidence"]
                )
                + "</h2>",
                return_val=False,
            )

            for reaction in collect_reactions_from_retrosynthesis(tree):
                if num_results < 4 or GLOBAL_SETTINGS["VERBOSE"] == False:
                    results[str(index)]["reactions"].append(reactions_text[i])
                output_text("<green> Reaction: </green>" + reactions_text[i], return_val=False)
                i = i + 1
                if GLOBAL_SETTINGS["display"] == "notebook":
                    display(Chem.Draw.ReactionToImage(reaction))
                else:
                    output_text("", return_val=False)

        save_result(
            create_analysis_record(product_smiles, "RXN", "Predict_Retrosynthesis", result_parameters, results),
            cmd_pointer=cmd_pointer,
        )
    except Exception as e:  # pylint: disable=broad-exception
        output_error(
            "The following error message was received while trying to display results: " + str(e),
            return_val=False,
        )
        return False
    i = 0

    if GLOBAL_SETTINGS["VERBOSE"] == False:
        return results
    else:
        return True
