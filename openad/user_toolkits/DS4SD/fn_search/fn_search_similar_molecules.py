"""searches for similar molecules"""
import numpy as np
from rdkit.Chem import PandasTools
from rdkit import Chem
import pandas as pd
from deepsearch.chemistry.queries.molecules import MoleculeQuery
from deepsearch.chemistry.queries.molecules import MolQueryType
from openad.helpers.output import output_text, output_success, output_error, output_table
from openad.helpers.output_msgs import msg
from openad.molecules.molecule_cache import create_analysis_record, save_result
from openad.molecules.mol_functions import canonical_smiles, valid_smiles
from openad.molecules.mol_commands import property_retrieve
from openad.app.global_var_lib import GLOBAL_SETTINGS

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
    except Exception as err:  # pylint: disable=broad-exception-caught
        output_error(msg("err_deepsearch", err), return_val=False)
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
        f"<h2>We found {len(results_table)} molecules similar to '<reset>{inputs['smiles']}</reset>'</h2>",
        return_val=False,
        pad_top=1,
    )
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
    if GLOBAL_SETTINGS["display"] == "notebook":
        from IPython.display import display

        if valid_smiles(inputs["smiles"]):
            try:
                smiles_mol = Chem.MolFromSmiles(inputs["smiles"])
            except Exception as err:  # pylint: disable= broad-exception-caught
                output_error("Error with rdkit verification of smiles:" + str(err), return_val=False)
                return False

            mol_img = Chem.Draw.MolToImage(smiles_mol, size=(200, 200))
            display(mol_img)

        df = pd.DataFrame(results_table)
        PandasTools.AddMoleculeColumnToFrame(df, smilesCol="SMILES")
        col = df.pop("ROMol")
        df.insert(0, col.name, col)
        col = df.pop("SMILES")
        df.insert(1, col.name, col)
        return df
    else:
        table = pd.DataFrame(results_table)
        output_table(table, tablefmt=_tableformat)
