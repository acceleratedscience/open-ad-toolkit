"""search for molecules in patents from list"""
import os
import numpy as np
from rdkit.Chem import PandasTools
from rdkit import Chem
import pandas as pd
from deepsearch.chemistry.queries.molecules import MoleculesInPatentsQuery
from openad.helpers.output import output_text, output_error, output_table, output_success
from openad.helpers.output_msgs import msg
from openad.app.global_var_lib import GLOBAL_SETTINGS


_tableformat = "simple"

""" needs to be migrated into Helper"""


def get_column_as_list_from_dataframe(a_dataframe, column_name) -> list:
    """returns a given column as a list object"""
    if column_name in a_dataframe:
        return a_dataframe[column_name].tolist()
    return []


def get_dataframe_from_file(cmd_pointer, filename):
    """returns a data frame from a csv file"""
    if not os.path.isfile(cmd_pointer.workspace_path(cmd_pointer.settings["workspace"].upper()) + "/" + filename):
        raise Exception("file does not exist")  # pylint: disable=broad-exception-raised
    else:
        df = pd.read_csv(cmd_pointer.workspace_path(cmd_pointer.settings["workspace"].upper()) + "/" + filename)
    return df


def search_molecules_in_patents(inputs: dict, cmd_pointer):
    """search patents for all mentioned molecules
    inputs: parser inputs from pyparsing
    cmd_pointer: pointer to runtime
    """
    api = cmd_pointer.login_settings["toolkits_api"][cmd_pointer.login_settings["toolkits"].index("DS4SD")]
    if isinstance(inputs["from_source"], dict) and inputs["from_source"]["from_list"] is not None:
        try:
            from_list = inputs["from_source"]["from_list"]
        except Exception:  # pylint: disable=broad-except
            output_error(
                "Unexpected pyparsing error. Please screenshot and report circumstance to OpenAD team",
                return_val=False,
            )
            output_error("Restart Notebook Kernel or application to proceed", return_val=False)
            return False

    elif "from_list" in inputs["from_source"][0]:
        try:
            from_list = inputs["from_source"][0]["from_list"]
        except Exception:  # pylint: disable=broad-except
            output_error(
                "Unexpected pyparsing error. Please screenshot and report circumstance to OpenAD team",
                return_val=False,
            )
            output_error("Restart Notebook Kernel or application to proceed", return_val=False)
            return False
    elif "from_dataframe" in inputs:
        try:
            react_frame = cmd_pointer.api_variables[inputs["from_dataframe"]]
            from_list = get_column_as_list_from_dataframe(react_frame, "PATENT ID")
            if from_list == []:
                from_list = get_column_as_list_from_dataframe(react_frame, "patent id")
            if from_list == []:
                raise Exception("No Patent ID found")  # pylint: disable=broad-exception-raised
        except Exception as err:  # pylint: disable=broad-except
            output_error(
                "Could not load valid list from dataframe column column 'PATENT ID' or 'patent id' " + str(err),
                return_val=False,
            )
            return True
    elif "from_file" in inputs:
        from_file = inputs["from_file"]
        try:
            react_frame = get_dataframe_from_file(cmd_pointer, from_file)
            from_list = get_column_as_list_from_dataframe(react_frame, "PATENT ID")
            if from_list == []:
                from_list = get_column_as_list_from_dataframe(react_frame, "patent id")
            if from_list == []:
                raise Exception("No Patent ID found")  # pylint: disable=broad-exception-raised
        except Exception:  # pylint: disable=broad-except
            output_error(
                "Could not load valid list from file column 'PATENT ID' or 'patent id' ",
                return_val=False,
            )
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

    except Exception as err:  # pylint: disable=broad-except
        output_error(["There was an error calling DeepSearch", err], return_val=False)
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
        output_success(msg("success_file_saved"), return_val=False)
    output_text(
        f"<h2>We found {len(results_table)} molecules that are mentioned in the following patents:</h2>",
        return_val=False,
        pad_top=1,
    )
    output_text("\n".join(from_list), return_val=False)

    if GLOBAL_SETTINGS["display"] == "notebook":
        df = pd.DataFrame(results_table)
        if len(df) > 0:
            PandasTools.AddMoleculeColumnToFrame(df, smilesCol="SMILES")
            col = df.pop("ROMol")
            df.insert(0, col.name, col)
            col = df.pop("SMILES")
            df.insert(1, col.name, col)
        return df
    else:
        table = pd.DataFrame(results_table)
        output_table(table, tablefmt=_tableformat)
