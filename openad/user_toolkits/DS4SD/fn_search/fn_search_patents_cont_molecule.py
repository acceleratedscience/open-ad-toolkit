# Example command:
# search for patents containing molecule 'CC(C)(c1ccccn1)C(CC(=O)O)Nc1nc(-c2c[nH]c3ncc(Cl)cc23)c(C#N)cc1F'
# search for patents containing molecule 'InChI=1S/C11H16N2OS/c1-4-8(14)7-9(12)11(2,3)10-13-5-6-15-10/h5-7H,4,12H2,1-3H3'

import numpy as np
import pandas as pd
from deepsearch.chemistry.queries.molecules import PatentsWithMoleculesQuery
from deepsearch.chemistry.queries.molecules import MolId, MolIdType
from openad.helpers.output import output_text, output_success, output_warning, output_error, output_table
from openad.helpers.output_msgs import msg
from openad.molecules.molecule_cache import create_analysis_record, save_result
from openad.molecules.mol_functions import canonical_smiles, valid_smiles, valid_inchi
from openad.molecules.mol_commands import property_retrieve

import os
import sys

parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, parent_dir)
from msgs import ds4sd_msg


def search_patents_cont_molecule(inputs: dict, cmd_pointer):
    """
    Searches patents that contain mentions of a given molecule.

    Parameters
    ----------
    inputs:
        Parser inputs from pyparsing.
    cmd_pointer:
        Pointer to runtime.
    """

    api = cmd_pointer.login_settings["toolkits_api"][cmd_pointer.login_settings["toolkits"].index("DS4SD")]
    result_type = ""
    try:
        if valid_smiles(inputs["smiles"]) is True:
            query = PatentsWithMoleculesQuery(
                molecules=[MolId(type=MolIdType.SMILES, value=canonical_smiles(inputs["smiles"]))],
                num_items=20,
            )

            result_type = "SMILES"
        elif valid_inchi(inputs["smiles"]) is True:
            query = PatentsWithMoleculesQuery(
                molecules=[MolId(type=MolIdType.INCHI, value=inputs["smiles"])],
                num_items=20,
            )
            result_type = "InChI"
        else:
            query = PatentsWithMoleculesQuery(
                molecules=[MolId(type=MolIdType.INCHIKEY, value=inputs["smiles"])],
                num_items=20,
            )
            result_type = "InChIKey"
            output_warning("String is not a SMILES or InChI, attempting an InChIKey search:", return_val=False)

        resp = api.queries.run(query)
        # raise Exception('This is a test error')
    except Exception as err:  # pylint: disable=broad-exception-caught
        output_error(ds4sd_msg("err_deepsearch", err), return_val=False)
        return False
    results_table = []

    for doc in resp.outputs["patents"]:
        result = {
            "PATENT ID": "",
        }
        for ident in doc["identifiers"]:
            if ident["type"] == "patentid":
                result["PATENT ID"] = ident["value"]
        results_table.append(result)

    if "save_as" in inputs:
        results_file = str(inputs["results_file"])
        df = pd.DataFrame(results_table)
        if not results_file.endswith(".csv"):
            results_file = results_file + ".csv"
        df.to_csv(
            cmd_pointer.workspace_path(cmd_pointer.settings["workspace"].upper()) + "/" + results_file, index=False
        )
        df = df.replace(np.nan, "", regex=True)
        output_success(msg("success_file_saved", results_file), return_val=False, pad_top=1, pad_btm=0)
    output_text(
        f"<bold>We found {len(results_table)} patents containing the requested {result_type}</bold>",
        return_val=False,
        pad_top=1,
    )
    output_text(inputs["smiles"], return_val=False)
    save_result(
        create_analysis_record(
            property_retrieve(inputs["smiles"], "canonical_smiles", cmd_pointer),
            "DS4SD",
            "Patents_For_Molecule",
            "",
            results_table,
        ),
        cmd_pointer=cmd_pointer,
    )

    return output_table(pd.DataFrame(results_table))
