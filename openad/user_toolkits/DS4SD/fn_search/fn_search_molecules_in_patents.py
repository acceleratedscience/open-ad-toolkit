# Example command:
# search for molecules in patents from list ['CN108473493B','US20190023713A1']

import os
import numpy as np
from rdkit.Chem import PandasTools
import pandas as pd
from deepsearch.chemistry.queries.molecules import MoleculesInPatentsQuery
from openad.app.global_var_lib import GLOBAL_SETTINGS
from openad.helpers.output import output_text, output_error, output_table, output_success
from openad.helpers.output_msgs import msg
from openad.helpers.h_data import col_from_df, csv_to_df
from openad.helpers.general import load_tk_module


def search_molecules_in_patents(inputs: dict, cmd_pointer):
    """
    Search for mentions of a given molecules in a list of patents.

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
    if isinstance(inputs["from_source"], dict) and inputs["from_source"]["from_list"] is not None:
        try:
            from_list = inputs["from_source"]["from_list"]
            # raise Exception('This is a test error')
        except Exception:  # pylint: disable=broad-except
            output_error(msg("err_pyparsing"), return_val=False)
            return False

    elif "from_list" in inputs["from_source"][0]:
        try:
            from_list = inputs["from_source"][0]["from_list"]
            # raise Exception('This is a test error')
        except Exception:  # pylint: disable=broad-except
            output_error(msg("err_pyparsing"), return_val=False)
            return False
    elif "from_dataframe" in inputs:
        try:
            react_frame = cmd_pointer.api_variables[inputs["from_dataframe"]]

            from_list = col_from_df(react_frame, "PATENT ID")
            if from_list == []:
                from_list = col_from_df(react_frame, "patent id")
            if from_list == []:
                raise Exception("No patent ID found")  # pylint: disable=broad-exception-raised
            # raise Exception('This is a test error')
        except Exception as err:  # pylint: disable=broad-except
            output_error(ds4sd_msg("err_patent_id", err), return_val=False)
            return True
    elif "from_file" in inputs:
        from_file = inputs["from_file"]
        try:
            react_frame = csv_to_df(cmd_pointer, from_file)
            from_list = col_from_df(react_frame, "PATENT ID")
            if from_list == []:
                from_list = col_from_df(react_frame, "patent id")
            if from_list == []:
                raise Exception("No patent ID found")  # pylint: disable=broad-exception-raised
            # raise Exception('This is a test error')
        except Exception as err:  # pylint: disable=broad-except
            output_error(ds4sd_msg("err_patent_id", err), return_val=False)
            return True
    try:
        query = MoleculesInPatentsQuery(
            patents=from_list,
            num_items=20,
        )

        resp = api.queries.run(query)
        results_table = []
        for row in resp.outputs["molecules"]:
            result = {
                "Id": row["persistent_id"],
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
        # raise Exception('This is a test error')
    except Exception as err:  # pylint: disable=broad-except
        output_error(ds4sd_msg("err_deepsearch", err), return_val=False)
        return False
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
        f"<bold>We found {len(results_table)} molecules that are mentioned in the following patents:</bold>",
        return_val=False,
        pad_top=1,
    )
    output_text("\n".join(from_list), return_val=False)

    df = pd.DataFrame(results_table)

    if GLOBAL_SETTINGS["display"] == "notebook":
        if len(df) > 0:
            PandasTools.AddMoleculeColumnToFrame(df, smilesCol="SMILES")
            col = df.pop("ROMol")
            df.insert(0, col.name, col)
            col = df.pop("SMILES")
            df.insert(1, col.name, col)
        # return output_table(df, is_data=True).data
        return df

    return output_table(df)
