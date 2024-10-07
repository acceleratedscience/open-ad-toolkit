"""
Functions to translate between different molecule and molecule set formats.
"""

import ast
import json
import pandas as pd
from rdkit import Chem, RDLogger

from openad.helpers.files import open_file
from openad.helpers.output import output_error
import openad.smols.smol_functions as smol_functions

# Suppress RDKit errors
RDLogger.DisableLog("rdApp.error")
RDLogger.DisableLog("rdApp.warning")


#
#


def smol2svg(inchi_or_smiles, highlight=None):
    """
    Takes an RDKit molecule object and returns an SVG string.

    Parameters
    ----------
    inchi_or_smiles: str
        Source option B: An InChI or SMILES identifier.
    highlight: str
        A SMARTS string to highlight a substructure in the molecule.
    """

    if not inchi_or_smiles:
        raise ValueError("Please provide inchi_or_smiles.")

    # Generate RDKit molecule object.
    mol_rdkit = Chem.MolFromInchi(inchi_or_smiles)
    if not mol_rdkit:
        mol_rdkit = Chem.MolFromSmiles(inchi_or_smiles)  # pylint: disable=no-member

    if highlight:
        substructure = Chem.MolFromSmarts(highlight)  # pylint: disable=no-member
        matches = mol_rdkit.GetSubstructMatches(substructure)

        # Flatten the tuple of tuples into a list of atom indices
        highlight_atoms = [atom_index for match in matches for atom_index in match]
    else:
        highlight_atoms = None

    mol_drawer = Chem.Draw.MolDraw2DSVG(400, 300)  # pylint: disable=no-member
    mol_drawer.DrawMolecule(mol_rdkit, highlightAtoms=highlight_atoms)
    mol_drawer.FinishDrawing()
    return mol_drawer.GetDrawingText()


def smol2mdl(smol=None, inchi_or_smiles=None, path=None):
    """
    Takes an RDKit molecule object OR an InChI/SMILES identifier and returns MDL data.

    An MDL molfile has the .mol extension and contains information about a molecule's atoms and bonds.
    It is used to store 2D and 3D structures of molecules, and required for our 3D visualization.
    https://en.wikipedia.org/wiki/Chemical_table_file#Molfile

    Note: an SDF file is a collection of MDL molfiles plus optional properties, separated by "$$$$"

    Parameters
    ----------
    smol: dict
        Source option A: An OpenAD small molecule dictionary (see OPENAD_SMOL_DICT).
    inchi_or_smiles: str
        Source option B: An InChI or SMILES identifier.
    path: str
        When path is defined, the MDL file will be written to disk at the specified location.
    """

    if not smol and not inchi_or_smiles:
        raise ValueError("Please provide smol or inchi_or_smiles.")

    inchi_or_smiles = smol["identifiers"]["inchi"] if smol else inchi_or_smiles

    # Generate RDKit molecule object.
    mol_rdkit = Chem.MolFromInchi(inchi_or_smiles)
    if not mol_rdkit:
        mol_rdkit = Chem.MolFromSmiles(inchi_or_smiles)  # pylint: disable=no-member

    # Add hydrogen atoms, which are displayed as spikes in the 3D viz.
    mol_rdkit = Chem.AddHs(mol_rdkit)  # pylint: disable=no-member

    # Generate 3D coordinates.
    Chem.rdDistGeom.EmbedMolecule(mol_rdkit)  # pylint: disable=no-member

    # Generate MDL data.
    mol_mdl = Chem.MolToMolBlock(mol_rdkit)  # pylint: disable=no-member

    # Write to disk
    if path:
        # Add smol properties
        if smol:
            for key, val in smol["properties"].items():
                val_print = "" if val is None else f"{val}"
                mol_rdkit.SetProp(key, val_print)
            mol_rdkit.SetProp("synonyms", f"{smol['synonyms']}")

        # # Print all properties - for debugging
        # props = mol_rdkit.GetPropsAsDict()
        # for key in props:
        #     print(">", key, props[key])

        # Write mdl to disk.
        with Chem.SDWriter(path) as writer:  # pylint: disable=no-member
            writer.write(mol_rdkit)

    # Return data
    else:
        # # We don't need the properties when returning data,
        # because this is only used for visualization. But for
        # the record, this is how we can add the properties:
        # if smol:
        #     props = smol["properties"].items()
        #     props_str = "\n" + "\n".join([f">  <{key}>\n{value}\n" for key, value in props])
        #     props_str = props_str + "\n" + f"<synonyms>\n{smol["synonyms"]}\n"
        #     mol_mdl += f"\n{props_str}\n$$$$\n"

        return mol_mdl


# Not used, for testing
def smol2xyz(mol_rdkit):
    """
    Takes an RDKit molecule object and returns it as XYZ data.
    """
    mol_xyz = Chem.rdmolfiles.MolToXYZBlock(mol_rdkit)
    return mol_xyz


# Not used, for testing
def smol2pdb(mol_rdkit):
    """
    Takes an RDKit molecule object and returns it as PDB data.
    """
    mol_pdb = Chem.rdmolfiles.MolToPDBBlock(mol_rdkit, flavor=32)
    return mol_pdb


# Not currently used.
# MolToMolBlock doesn't support parameters for the SDF file, so this code is not correct.
def dataframe2sdf(df):
    """
    Takes a dataframe with a SMILES or InChI column and returns SDF data.
    The other columns will be included as properties.

    Parameters
    ----------
    df: DataFrame
        A pandas DataFrame with a SMILES or InChI column.
    """

    # This allows us to do a case-insensitive scan for the InChI or SMILES column.
    cols_lowercase = [col.lower() for col in df.columns]

    if "inchi" in cols_lowercase:
        index = cols_lowercase.index("inchi")
        key = df.columns[index]
        key_type = "inchi"
    elif "smiles" in cols_lowercase:
        index = cols_lowercase.index("smiles")
        key = df.columns[index]
        key_type = "smiles"
    else:
        return None

    # Convert the molecules to SDF format
    sdf_data = ""
    for i, row in df.iterrows():
        if key_type == "inchi":
            mol_rdkit = Chem.MolFromInchi(row[key])  # pylint: disable=no-member
        elif key_type == "smiles":
            mol_rdkit = Chem.MolFromSmiles(row[key])  # pylint: disable=no-member

        if mol_rdkit is not None:
            # Add all other dataframe columns as properties to the SDF data,
            # unless they're other identifiers, in which case they will be ignored.
            for col in df.columns:
                # if col.lower() not in ["inchi", "smiles"]: # %%
                mol_rdkit.SetProp(col, str(row[col]))

            sdf_data += Chem.MolToMolBlock(mol_rdkit) + "\n$$$$\n"  # pylint: disable=no-member

    return sdf_data


def dataframe2molset(df):
    """
    Takes a dataframe with a SMILES or InChI column and returns a molset dictionary.
    The other columns will be included as properties.

    Parameters
    ----------
    df: DataFrame
        A pandas DataFrame with a SMILES or InChI column.
    """

    # This allows us to do a case-insensitive scan for the InChI or SMILES column.
    cols_lowercase = [col.lower() for col in df.columns]

    if "inchi" in cols_lowercase:
        identifier_index = cols_lowercase.index("inchi")
        identifier = df.columns[identifier_index]
    elif "smiles" in cols_lowercase:
        identifier_index = cols_lowercase.index("smiles")
        identifier = df.columns[identifier_index]
    elif "canonical_smiles" in cols_lowercase:
        identifier_index = cols_lowercase.index("canonical_smiles")
        identifier = df.columns[identifier_index]
    elif "isomeric_smiles" in cols_lowercase:
        identifier_index = cols_lowercase.index("isomeric_smiles")
        identifier = df.columns[identifier_index]
    else:
        return None

    # Convert the molecules to SDF format
    molset = []
    for i, row in df.iterrows():
        # Get name field regardless of case.
        name = next((value for key, value in row.items() if key.lower() == "name"), None)

        # Create molecule
        mol_dict = smol_functions.new_smol(row[identifier], name=name)

        if mol_dict is not None:
            # Add all other dataframe columns as properties to the SDF data,
            # unless they're other identifiers, in which case they will be ignored. %%
            for col in df.columns:
                mol_dict["properties"][col] = str(row[col])

            # Add index
            mol_dict["index"] = i + 1

            molset.append(mol_dict)

    return molset


def molset2dataframe(molset, remove_invalid_mols=False, include_romol=False):
    """
    Takes a molset dictionary and returns an Pandas dataframe.

    Parameters
    ----------
    molset: list
        A list of molecule objects, our OpenAD molset format.
    remove_invalid_mols: bool
        Unless set to True, the function will fail and return a list
        of invalid molecules if any of the mols in the molset cannot
        be parsed by RDKit. This information is used in the GUI to then
        display the list of molecules that will be removed if the user
        chooses to proceed. After confirming, the function will run again
        with remove_invalid_molsset to True.
    include_romol: bool
        If set to True, the RDKit molecule object will be included in the
        dataframe. We need it to transform the dataframe to an SDF of MDL file.
    """

    # Flatten the molset into a list of dictionaries.
    data = []
    invalid = []
    for i, smol in enumerate(molset):
        # Create RDKit molecule object (ROMol)
        if smol["identifiers"].get("inchi"):
            mol_rdkit = Chem.MolFromInchi(smol["identifiers"]["inchi"])
        else:
            smiles = smol_functions.get_best_available_smiles(smol)
            if smiles:
                mol_rdkit = Chem.MolFromSmiles(smiles)  # pylint: disable=no-member

        # Store mols that failed to parse.
        if not mol_rdkit:
            output_error("Failed to parse:", i, smiles)  # Keep this
            invalid.append(i)
            continue

        # Add ROMol the row.
        if include_romol:
            row = {"ROMol": mol_rdkit}
        else:
            row = {}

        # Add identifiers to the row.
        for key in smol["identifiers"]:
            if smol["identifiers"][key]:
                row[key] = smol["identifiers"][key]

        # Add properties to the row.
        for key in smol["properties"]:
            if smol["properties"][key] is None:  # To avoid None values to be stored as "None" string
                row[key] = ""
            else:
                row[key] = smol["properties"][key]

        # Add synonyms to the row.
        if smol["synonyms"] and len(smol["synonyms"]) > 0:
            row["synonyms"] = "\n".join(smol["synonyms"])
        data.append(row)

    # Return list of failures unless user has
    # explicitly requested to remove them.
    if invalid and not remove_invalid_mols:
        invalid_mols = [molset[i] for i in invalid]
        raise ValueError(f"Failed to parse {len(invalid)} molecules", invalid_mols)

    # Turn data into dataframe
    df = pd.DataFrame(data)
    df = df.fillna("")  # Replace NaN with empty string
    return df


def write_dataframe2sdf(df, destination_path):
    """
    Takes a dataframe with a molecule column and writes it to an SDF file.

    Parameters
    ----------
    df: DataFrame
        A pandas DataFrame with a molecule column.
    destination_path: str
        The path to save the SDF file to.
    """

    if "ROMol" not in df.columns:
        raise ValueError("Dataframe does not contain a 'ROMol' column")

    try:
        Chem.PandasTools.WriteSDF(df, destination_path, molColName="ROMol", properties=list(df.columns), idName="RowID")
    except Exception as err:
        raise RuntimeError(f"Failed to write dataframe to SDF: {err}")


def write_dataframe2csv(df, destination_path):
    """
    Takes any dataframe and writes it to a CSV file.

    Parameters
    ----------
    df: DataFrame
        A pandas DataFrame.
    destination_path: str
        The path to save the CSV file to.
    """

    # Remove the ROMol (RDKit molecule object) column if it exists.
    if "ROMol" in df.columns:
        df = df.drop(columns=["ROMol"])

    try:
        df.to_csv(destination_path, index=False)
    except Exception as err:
        raise RuntimeError(f"Failed to write dataframe to CSV: {err}")


# Return molset from SMILES file
def smiles_path2molset(path_absolute):
    """
    Takes the content of a .smi file and returns a molset dictionary.
    Specs for .smi files: http://opensmiles.org/opensmiles.html - 4.5

    This takes about 3 seconds per 10,000 molecules on an Apple M2 with 16GB or memory.
    """

    # Read file's content
    data, err_code = open_file(path_absolute, return_err="code")
    if err_code:
        return None, err_code

    # Parse SMILES
    smiles_list = data.splitlines()
    # Ignore any properties that may be listed after the SMILES string.
    smiles_list = [smiles.split(" ")[0] for smiles in smiles_list if smiles]
    molset = []
    for i, smiles in enumerate(smiles_list):
        mol = smol_functions.new_smol(smiles)
        if not mol:
            mol = {
                "identifiers": {"canonical_smiles": smiles},
                "properties": {},
            }

        mol["index"] = i + 1
        molset.append(mol)

    return molset, None


def sdf_path2molset(sdf_path):
    """
    Takes the content of an .sdf file and returns a molset dictionary.

    Used to open SDF files.
    """

    # This lets us parse all the properties back to their original types,
    # since SDF stores them all as string. However this is can cause
    # unexpected results, without providing any real benefit. For example,
    # The cactvs_fingerprint property stores a long binary string, which
    # gets converted into a number, which then becomes "Infinite" after parsing.
    def _try_parse_json(value):
        try:
            return json.loads(value)
        except json.JSONDecodeError:
            try:
                return ast.literal_eval(value)
            except (ValueError, SyntaxError):
                return value

    try:
        mols_rdkit = Chem.SDMolSupplier(sdf_path)  # pylint: disable=no-member
        molset = []
        for i, mol_rdkit in enumerate(mols_rdkit):
            mol_dict = smol_functions.new_smol(mol_rdkit=mol_rdkit)
            mol_dict["index"] = i + 1
            molset.append(mol_dict)
        return molset, None
    except Exception as err:
        return None, err


def csv_path2molset(csv_path):
    """
    Takes the content of a .csv file and returns a molset dictionary.

    Used to open CSV files.
    """

    # Read CSV data
    try:
        df = pd.read_csv(csv_path)
        # Convert to molset
        molset = dataframe2molset(df)
        return molset, None
    except Exception as err:
        return None, err


def mdl_path2smol(mdl_path):
    """
    Takes the content of a .mol file and returns a small molecule dictionary.

    Used to open MDL (.mol) files.
    """

    # Read MDL data
    supplier = Chem.SDMolSupplier(mdl_path)  # pylint: disable=no-member
    mol_rdkit = next(supplier)

    if mol_rdkit is None:
        return None, "unknown"

    # Translate into OpenAD smol dict
    mol_dict = smol_functions.new_smol(mol_rdkit=mol_rdkit)
    return mol_dict, None
