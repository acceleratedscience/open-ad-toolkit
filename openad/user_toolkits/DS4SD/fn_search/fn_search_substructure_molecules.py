import numpy as np
import pandas as pd
from rdkit.Chem import PandasTools
from rdkit import Chem
from deepsearch.chemistry.queries.molecules import MoleculeQuery
from deepsearch.chemistry.queries.molecules import MolQueryType
from openad.helpers.output import output_text, output_success, output_error, output_table
from openad.helpers.output_msgs import msg
from openad.app.global_var_lib import GLOBAL_SETTINGS

_tableformat = "simple"


# needs to be migrated into Helper
def valid_smiles(input_molecule) -> bool:
    """tests to see if a input molecule is valid smiles definition
    input_molecule: smiles string"""
    from rdkit import rdBase

    blocker = rdBase.BlockLogs()  # pylint: disable=c-extension-no-member
    try:
        m = Chem.MolFromSmiles(input_molecule, sanitize=False)  # pylint: disable=no-member
    except:
        return False
    if m is None:
        return False
    else:
        try:
            Chem.SanitizeMol(m)  # pylint: disable=no-member
        except Exception:  # pylint: disable=broad-exception-caught
            return False
    return True


def search_substructure_molecules(inputs: dict, cmd_pointer):
    api = cmd_pointer.login_settings["toolkits_api"][cmd_pointer.login_settings["toolkits"].index("DS4SD")]
    try:
        query = MoleculeQuery(
            query=inputs["smiles"],
            query_type=MolQueryType.SUBSTRUCTURE,
        )
        resp = api.queries.run(query)
    except Exception as e:  # pylint: disable=broad-exception-caught
        output_error("Error in calling deepsearch:" + str(e), return_val=False)
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
        output_success(msg("success_file_saved"), return_val=False)
    output_text(
        "\n<h2>Search Results for smiles containing substructures of : </h2> ",
        return_val=False,
    )
    output_text(inputs["smiles"], return_val=False)

    if GLOBAL_SETTINGS["display"] == "notebook":
        from IPython.display import display

        if valid_smiles(inputs["smiles"]):
            try:
                smiles_mol = Chem.MolFromSmiles(inputs["smiles"])
            except Exception as e:  # pylint: disable=broad-exception-caught
                output_error("Error with rdkit verification of smiles:" + str(e), return_val=False)
                return False

            display(smiles_mol)
        df = pd.DataFrame(results_table)
        col = df.pop("SMILES")
        df.insert(0, col.name, col)
        PandasTools.AddMoleculeColumnToFrame(df, smilesCol="SMILES")
        col = df.pop("ROMol")
        df.insert(1, col.name, col)
        return df
    else:
        table = pd.DataFrame(results_table)
        output_table(table, tablefmt=_tableformat)
