""" Searches Patents for occurences of a molecule"""
import numpy as np
import pandas as pd
from deepsearch.chemistry.queries.molecules import PatentsWithMoleculesQuery
from deepsearch.chemistry.queries.molecules import MolId, MolIdType
from rdkit import Chem
from openad.helpers.output import output_table
from openad.helpers.output import output_error
from openad.helpers.output import output_text
from openad.helpers.output import output_warning

_tableformat = "simple"


# needs to be migrated into Helper
def valid_smiles(input_molecule) -> bool:
    """tests to see if a input molecule is valid smiles definition
    input_molecule: smiles string"""
    from rdkit import rdBase

    blocker = rdBase.BlockLogs()  # pylint: disable=c-extension-no-member
    m = Chem.MolFromSmiles(input_molecule, sanitize=False)  # pylint: disable=no-member
    if m is None:
        return False
    else:
        try:
            Chem.SanitizeMol(m)  # pylint: disable=no-member
        except Exception:  # pylint: disable=broad-exception-caught
            return False
    return True


def valid_inchi(input_molecule) -> bool:
    """tests to see if a input molecule is valid inchi definition
    input_molecule: inchi string"""
    from rdkit import rdBase

    blocker = rdBase.BlockLogs()  # pylint: disable=c-extension-no-member
    m = Chem.inchi.InchiToInchiKey(input_molecule)
    if m is None:
        return False
    else:
        return True


def search_patent_molecule(inputs: dict, cmd_pointer):
    """Searches Patents for occurences of a molecule
    inputs: parser inputs from pyparsing
    cmd_pointer: pointer to runtime
    """

    api = cmd_pointer.login_settings["toolkits_api"][cmd_pointer.login_settings["toolkits"].index("DS4SD")]
    result_type = ""
    try:
        if valid_smiles(inputs["smiles"]) is True:
            query = PatentsWithMoleculesQuery(
                molecules=[MolId(type=MolIdType.SMILES, value=inputs["smiles"])],
                num_items=20,
            )
            result_type = "SMILES"
        elif valid_inchi(inputs["smiles"]) is True:
            query = PatentsWithMoleculesQuery(
                molecules=[MolId(type=MolIdType.INCHI, value=inputs["smiles"])],
                num_items=20,
            )
            result_type = "Inchi"
        else:
            query = PatentsWithMoleculesQuery(
                molecules=[MolId(type=MolIdType.INCHIKEY, value=inputs["smiles"])],
                num_items=20,
            )
            result_type = "Inchi Key"
            output_warning(
                "String is Not a Smiles or Inchi attemnting Inchikey search: ",
                cmd_pointer=cmd_pointer,
                return_val=False,
            )
        resp = api.queries.run(query)
    except Exception as e:  # pylint: disable=broad-exception-caught
        output_error("Error in calling deepsearch:" + str(e), cmd_pointer=cmd_pointer, return_val=False)
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
        output_text(
            "\n <success>File successfully saved to workspace.</success>", cmd_pointer=cmd_pointer, return_val=False
        )
    output_text("", cmd_pointer=cmd_pointer, return_val=False)
    output_text(
        " <h2>  Patent Search Results for " + result_type + " molecule:</h2> ",
        cmd_pointer=cmd_pointer,
        return_val=False,
    )
    output_text(inputs["smiles"], cmd_pointer=cmd_pointer, return_val=False)

    table = pd.DataFrame(results_table)
    return output_table(table, cmd_pointer, tablefmt=_tableformat)
