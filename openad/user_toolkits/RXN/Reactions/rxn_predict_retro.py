"""This library is for implementing the predict retrosynthesis function from RXN"""

from typing import Dict, List
import importlib.util as ilu
from openad.helpers.output import output_text
from openad.helpers.output import output_error
from rdkit import Chem
from rdkit.Chem import AllChem
from time import sleep


_tableformat = "simple"


def get_include_lib(cmd_pointer):
    """load the rxn include library functions"""
    folder = cmd_pointer.toolkit_dir + "/RXN" + "/rxn_include.py"
    file = "rxn_include"
    spec = ilu.spec_from_file_location(file, folder)
    rxn = ilu.module_from_spec(spec)
    spec.loader.exec_module(rxn)
    rxn_helper = rxn.rxn_helper()
    return rxn_helper


def get_reaction_from_smiles(reaction_smiles: str) -> Chem.rdChemReactions.ChemicalReaction:
    """get reaction image"""
    return AllChem.ReactionFromSmarts(reaction_smiles, useSmiles=True)  # pylint: disable=no-member


def collect_reactions_from_retrosynthesis(tree: Dict) -> List[str]:
    """collect all reactions from retrosynthesis tree"""
    reactions = []
    if "children" in tree and len(tree["children"]):
        reactions.append(
            AllChem.ReactionFromSmarts(
                "{}>>{}".format(".".join([node["smiles"] for node in tree["children"]]), tree["smiles"]), useSmiles=True
            )
        )  # pylint: disable=no-member
    for node in tree["children"]:
        reactions.extend(collect_reactions_from_retrosynthesis(node))
    return reactions


def collect_reactions_from_retrosynthesis_text(tree: Dict) -> List[str]:
    """collect all reactions from retrosynthesis tree"""
    reactions = []
    if "children" in tree and len(tree["children"]):
        reactions.append(
            "{} --->> {}".format(" + ".join([node["smiles"] for node in tree["children"]]), tree["smiles"])
        )

    for node in tree["children"]:
        reactions.extend(collect_reactions_from_retrosynthesis_text(node))

    return reactions


def predict_retro(inputs: dict, cmd_pointer):
    """Perform RXN Predict Retro Synthesis"""
    rxn_helper = get_include_lib(cmd_pointer)
    rxn_helper.sync_up_workspace_name(cmd_pointer)
    rxn_helper.get_current_project(cmd_pointer)
    if cmd_pointer.notebook_mode is True:
        from IPython.display import display  # pylint: disable=import-outside-toplevel
        from halo import HaloNotebook as Halo  # pylint: disable=import-outside-toplevel
    else:
        from halo import Halo  # pylint: disable=import-outside-toplevel

    class Spinner(Halo):
        "contextual spinner"

        def __init__(self):
            # Alternative spinners:
            # simpleDotsScrolling, interval=100
            super().__init__(spinner="dots", color="white")

    val = "val"
    availability_pricing_threshold = 0
    available_smiles = None
    exclude_smiles = None
    exclude_substructures = None
    exclude_target_molecule = False
    fap = 0.6
    max_steps = 3
    nbeams = 10
    pruning_steps = 2
    ai_model = "2020-07-01"

    product_smiles = inputs["molecule"]
    #######################

    if not rxn_helper.valid_smiles(str(product_smiles)):
        output_error(" Invalid Smiles Supplied.", cmd_pointer=cmd_pointer, return_val=False)
        return False

    if cmd_pointer.notebook_mode is True:
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
        output_text("<success>Target Molecule:</success> " + product_smiles, cmd_pointer=cmd_pointer, return_val=False)
    #######################

    if "availability_pricing_threshold" in inputs:
        availability_pricing_threshold = int(inputs["availability_pricing_threshold"][val])
    if "available_smiles" in inputs:
        available_smiles = inputs["available_smiles"][val]
    if "exclude_smiles" in inputs:
        exclude_smiles = inputs["exclude_smiles"][val]
    if "exclude_substructures" in inputs:
        exclude_substructures = inputs["exclude_substructures"][val]
    if "exclude_target_molecule" in inputs:
        if inputs["exclude_substructures"][val].upper() == "TRUE":
            exclude_target_molecule = True
    if "fap" in inputs:
        fap = float(inputs["fap"][val])
    if "max_steps" in inputs:
        max_steps = int(inputs["max_steps"][val])
    if "nbeams" in inputs:
        nbeams = int(inputs["nbeams"][val])
    if "pruning_steps" in inputs:
        pruning_steps = int(inputs["pruning_steps"][val])
    if "ai_model" in inputs:
        ai_model = inputs["ai_model"][val]
    rxn4chemistry_wrapper = cmd_pointer.login_settings["client"][cmd_pointer.login_settings["toolkits"].index("RXN")]

    # Prepare the data query
    print("\n")
    newspin = Spinner()
    newspin.start("Starting Retrosynthesis")
    try:
        retries = 0
        status = False
        predict_retro_response = None
        while retries < 10 and status is False:
            try:
                newspin.text = "Submitting Retrosynthesis "
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

        while status != "SUCCESS":
            try:
                newspin.text = "Processing Retrosynthesis :" + status
                predict_automatic_retrosynthesis_results = (
                    rxn4chemistry_wrapper.get_predict_automatic_retrosynthesis_results(
                        predict_retro_response["prediction_id"]
                    )
                )
                if predict_retro_response["response"]["payload"] is None:
                    output_text(
                        "<h2>No Result:</h2>  Unable to find path for  " + product_smiles,
                        cmd_pointer=cmd_pointer,
                        return_val=False,
                    )
                    if cmd_pointer.notebook_mode is True:
                        return
                    else:
                        return False
                status = predict_automatic_retrosynthesis_results["status"]
                sleep(5)

            except Exception as e:  # pylint: disable=broad-exception-caught
                retries = retries + 1
                sleep(15)
                newspin.text = "Processing Retrosynthesis: Waiting"
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

    try:
        for index, tree in enumerate(predict_automatic_retrosynthesis_results["retrosynthetic_paths"]):
            for reaction in collect_reactions_from_retrosynthesis_text(tree):
                reactions_text.append(str(reaction))
    except Exception as e:  # pylint: disable=broad-exception-caught
        newspin.fail("Unable to Process")
        newspin.stop()
        raise Exception(
            "The following Error message was received while trying to process results:" + str(e)
        ) from e  # pylint: disable=broad-exception-raised
    i = 0
    try:
        newspin.succeed("Finished Processing")
        newspin.start()
        newspin.stop()
        for index, tree in enumerate(predict_automatic_retrosynthesis_results["retrosynthetic_paths"]):
            output_text(
                "",
                cmd_pointer=cmd_pointer,
                return_val=False,
            )
            output_text(
                "<h2> <success>Showing path </success> {} <success> with confidence </success>{}:".format(
                    index, tree["confidence"]
                )
                + "</h2>",
                cmd_pointer=cmd_pointer,
                return_val=False,
            )
            for reaction in collect_reactions_from_retrosynthesis(tree):
                output_text(
                    "<success> Reaction: </success>" + reactions_text[i], cmd_pointer=cmd_pointer, return_val=False
                )
                if cmd_pointer.notebook_mode is True:
                    display(Chem.Draw.ReactionToImage(reaction))
                else:
                    output_text("", cmd_pointer=cmd_pointer, return_val=False)

    except Exception as e:  # pylint: disable=broad-exception
        output_error(
            "The following error message was received while trying to display results: " + str(e),
            cmd_pointer=cmd_pointer,
            return_val=False,
        )
        return False
    i = 0
    return True
