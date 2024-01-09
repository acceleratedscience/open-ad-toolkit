""" Searches Patents for occurences of a molecule"""
import numpy as np
import pandas as pd
from deepsearch.chemistry.queries.molecules import PatentsWithMoleculesQuery
from deepsearch.chemistry.queries.molecules import MolId, MolIdType
from rdkit import Chem
from openad.helpers.output import output_text, output_success, output_warning, output_error, output_table
from openad.helpers.output_msgs import msg
from openad.molecules.molecule_cache import create_analysis_record, save_result
from openad.molecules.mol_functions import canonical_smiles, valid_smiles
from openad.molecules.mol_commands import property_retrieve

_tableformat = "simple"


# needs to be migrated into Helper


def valid_inchi(input_molecule) -> bool:
    """tests to see if a input molecule is valid inchi definition
    input_molecule: inchi string"""
    from rdkit import rdBase

    blocker = rdBase.BlockLogs()  # pylint: disable=c-extension-no-member
    try:
        m = Chem.inchi.InchiToInchiKey(input_molecule)
    except:
        return False
    if m is None:
        return False
    else:
        return True


def search_patents_cont_molecule(inputs: dict, cmd_pointer):
    """Searches Patents for occurences of a molecule
    inputs: parser inputs from pyparsing
    cmd_pointer: pointer to runtime
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
            output_warning("String is Not a Smiles or Inchi attemnting Inchikey search: ", return_val=False)

        resp = api.queries.run(query)

    except Exception as err:  # pylint: disable=broad-exception-caught
        output_error(["There was an error calling DeepSearch", err], return_val=False)
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
        output_success(msg("success_file_saved"), return_val=False)
    output_text(
        f"<h2>We found {len(results_table)} patents containing the requested {result_type}</h2>", return_val=False
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

    table = pd.DataFrame(results_table)
    return output_table(table, tablefmt=_tableformat)
