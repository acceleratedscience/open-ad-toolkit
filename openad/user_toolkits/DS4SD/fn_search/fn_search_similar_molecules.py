# Example command:
# search for similar molecules to 'C1(C(=C)C([O-])C1C)=O'

import numpy as np
from rdkit.Chem import PandasTools
from rdkit import Chem
import pandas as pd
from deepsearch.chemistry.queries.molecules import MoleculeQuery
from deepsearch.chemistry.queries.molecules import MolQueryType
from openad.smols.smol_cache import create_analysis_record, save_result
from openad.smols.smol_functions import valid_smiles
from openad.app.global_var_lib import GLOBAL_SETTINGS
from openad.helpers.output import output_text, output_success, output_error, output_table
from openad.helpers.output_msgs import msg
from openad.helpers.general import load_tk_module
from openad.smols.smol_functions import canonicalize


def search_similar_molecules(inputs: dict, cmd_pointer):
    """
    Search for molecules similar to a given molecule.

    Parameters
    ----------
    inputs:
        Parser inputs from pyparsing.
    cmd_pointer:
        Pointer to runtime.
    """

    # Load module from the toolkit folder.
    ds4sd_msg = load_tk_module(cmd_pointer, "DS4SD", "msgs", "ds4sd_msg")

    api = cmd_pointer.login_settings["toolkits_api"][cmd_pointer.login_settings["toolkits"].index("DS4SD")]
    try:
        query = MoleculeQuery(
            query=canonicalize(inputs["smiles"]),
            query_type=MolQueryType.SIMILARITY,
        )

        resp = api.queries.run(query)
        # raise Exception('This is a test error')
    except Exception as err:  # pylint: disable=broad-exception-caught
        output_error(ds4sd_msg("err_deepsearch", err), return_val=False)
        return False
    results_table = []
    for row in resp.outputs["molecules"]:
        result = {
            "id": row["persistent_id"],
            "SMILES": "",
            "InChIKey": "",
            "InChI": "",
        }
        for ref in row["identifiers"]:
            if ref["type"] == "smiles":
                result["SMILES"] = ref["value"]
            if ref["type"] == "inchikey":
                result["InChIKey"] = ref["value"]
            if ref["type"] == "inchi":
                result["InChI"] = ref["value"]
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
        output_success(msg("success_file_saved", results_file), return_val=False, pad_top=1)

    output_text(
        f"<bold>We found {len(results_table)} molecules similar to the provided SMILES</bold>",
        return_val=False,
        pad_top=1,
    )
    output_text(inputs["smiles"], return_val=False)
    save_result(
        create_analysis_record(
            canonicalize(inputs["smiles"]),
            "DS4SD",
            "Similar_Molecules",
            "",
            results_table,
        ),
        cmd_pointer=cmd_pointer,
    )

    df = pd.DataFrame(results_table)

    if GLOBAL_SETTINGS["display"] == "notebook":
        from IPython.display import display

        if valid_smiles(inputs["smiles"]):
            try:
                smiles_mol = Chem.MolFromSmiles(inputs["smiles"])
                # raise Exception('This is a test error')
            except Exception as err:  # pylint: disable= broad-exception-caught
                output_error(ds4sd_msg("err_rdkit_smiles", err), return_val=False)
                return False

            mol_img = Chem.Draw.MolToImage(smiles_mol, size=(200, 200))
            display(mol_img)

        PandasTools.AddMoleculeColumnToFrame(df, smilesCol="SMILES")
        col = df.pop("ROMol")
        df.insert(0, col.name, col)
        col = df.pop("SMILES")
        df.insert(1, col.name, col)
        # return output_table(df, is_data=True).data
        return df

    return output_table(df, is_data=True)
