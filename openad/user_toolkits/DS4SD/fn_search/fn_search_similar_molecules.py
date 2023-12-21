"""searches for similar molecules"""
import numpy as np
from rdkit.Chem import PandasTools
from rdkit import Chem
import pandas as pd
from deepsearch.chemistry.queries.molecules import MoleculeQuery
from deepsearch.chemistry.queries.molecules import MolQueryType
from openad.helpers.output import output_table
from openad.helpers.output import output_error
from openad.helpers.output import output_text
from openad.molecules.molecule_cache import create_analysis_record, save_result
from openad.molecules.mol_functions import canonical_smiles, valid_smiles
from openad.molecules.mol_commands import property_retrieve

_tableformat = "simple"


# needs to be migrated into Helper


def search_similar_molecules(inputs: dict, cmd_pointer):
    """search for molecules similar to a specified Molecule
    inputs: parser inputs from pyparsing
    cmd_pointer: pointer to runtime
    """
    api = cmd_pointer.login_settings["toolkits_api"][cmd_pointer.login_settings["toolkits"].index("DS4SD")]
    try:
        query = MoleculeQuery(
            query=inputs["smiles"],
            query_type=MolQueryType.SIMILARITY,
        )

        resp = api.queries.run(query)
    except Exception as e:  # pylint: disable=broad-exception-caught
        output_error("Error in calling deepsearch:" + str(e), cmd_pointer=cmd_pointer, return_val=False)
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
        output_text(
            " <green> File successfully saved to workspace. </green>", cmd_pointer=cmd_pointer, return_val=False
        )
    output_text(" ", cmd_pointer=cmd_pointer, return_val=False)
    output_text(
        "<h2>  Similarity Search Results for smiles molecule: </h2>  ", cmd_pointer=cmd_pointer, return_val=False
    )
    output_text(inputs["smiles"], cmd_pointer=cmd_pointer, return_val=False)
    save_result(
        create_analysis_record(
            inputs["smiles"],
            "DS4SD",
            "Similar_Molecules",
            "",
            results_table,
        ),
        cmd_pointer=cmd_pointer,
    )
    if cmd_pointer.notebook_mode is True:
        from IPython.display import display

        if valid_smiles(inputs["smiles"]):
            try:
                smiles_mol = Chem.MolFromSmiles(inputs["smiles"])
            except Exception as e:  # pylint: disable= broad-exception-caught
                output_error(
                    "Error with rdkit verification of smiles:" + str(e), cmd_pointer=cmd_pointer, return_val=False
                )
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
        output_table(table, cmd_pointer, tablefmt=_tableformat)
