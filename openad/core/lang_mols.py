"""Commands to launch molecule grid and molecule viewer."""

import pandas as pd
from rdkit import Chem, RDLogger
from openad.flask_apps import launcher
from openad.flask_apps.molsgrid.routes import fetchRoutesMolsGrid
from openad.flask_apps.molviewer.routes import fetchRoutesMolViewer
from openad.helpers.general import open_file
from openad.helpers.output import msg, output_error

# Disable RDKit warnings
logger = RDLogger.logger()
logger.setLevel(RDLogger.CRITICAL)


# Launch molecule grid.
def show_molsgrid(cmd_pointer, parser):
    # Load routes and launch browser UI.
    routes, the_mols2grid = fetchRoutesMolsGrid(cmd_pointer, parser)
    if cmd_pointer.notebook_mode:
        return the_mols2grid
    else:
        launcher.launch(cmd_pointer, routes, "molsgrid")


# Launch molecule viewer.
def show_mol(cmd_pointer, parser):
    if not "input_str" in parser.as_dict():
        output_error(msg("input_str_missing"), pad=1)
        return

    # Identify input type.
    input_str = parser.as_dict()["input_str"]

    # mol_dict = our own molecule JSON format.
    # mol_obj = RDKit molecule object.
    mol_dict, mol_obj = _parse_mol(cmd_pointer, input_str)

    # Abort if needed.
    if not mol_dict:
        return

    _complete_identifiers(mol_dict, mol_obj)  # Fill in missing identifiers: InChI, InChIKey, SMILES.
    svg = _renderSVG(mol_obj)
    mol_sdf = _mol2sdf(mol_obj)

    # Load routes and launch browser UI.
    routes = fetchRoutesMolViewer(mol_dict, mol_sdf, svg)

    if cmd_pointer.notebook_mode:
        # Jupyter
        launcher.launch(cmd_pointer, routes, "molviewer")
    else:
        # CLI
        launcher.launch(cmd_pointer, routes, "molviewer")


# Take the input string and parse it into a mol_obj and mol_dict,
# regardless if it's a path, InChI, InChIKey, or SMILES.
def _parse_mol(cmd_pointer, input_str):
    input_type = __identitify_input_type(input_str)
    mol_dict = None
    mol_obj = None

    if input_type == "InChI":
        mol_obj = Chem.MolFromInchi(input_str)
        mol_dict = {"InChI": input_str}
    # elif input_type == "InChIKey":
    #     # We could check pubchem for this compound...
    #     # https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/BPGDAMSIGCZZLK-UHFFFAOYSA-N/SDF
    elif input_type == "SMILES":
        mol_obj = Chem.MolFromSmiles(input_str)
        mol_dict = {"SMILES": input_str}
    elif input_type == "path":
        # Compile complete file path.
        workspace_path = cmd_pointer.workspace_path(cmd_pointer.settings["workspace"].upper()) + "/"
        file_path = workspace_path + "_mols/" + input_str

        # Turn supported file formats into a dictionary.
        mol_dict = __file2dict(file_path)

        # Use InChI or SMILES to create mol_obj.
        if mol_dict:
            mol_obj = __create_mol_obj(mol_dict)

            # Abort if needed.
            if not "SMILES" in mol_dict and not "InChI" in mol_dict:
                output_error(msg("identifier_missing", file_path, split=True), pad=1)
                mol_dict = None

    return mol_dict, mol_obj


# Figure out what type of input we're dealing with.
def __identitify_input_type(input_str):
    # Check if it's InChI
    if input_str.startswith("InChI="):
        try:
            _ = Chem.MolFromInchi(input_str)
            if _:
                return "InChI"
        except Exception:
            pass

    # # Check if it's InChIKey
    # if len(input_str) == 27 and "-" in input_str and input_str.replace("-", "").isalnum():
    #     return "InChIKey"

    # Check if it's SMILES
    try:
        _ = Chem.MolFromSmiles(input_str)
        if _:
            return "SMILES"
    except Exception:
        pass

    # Assume it's a path
    return "path"


# Read json/sdf file and return a dictionary of molecule properties.
def __file2dict(file_path):
    ext = file_path.split(".")[-1]
    mol_dict = None

    # JSON
    if ext == "json" or ext == "cjson":
        mol_dict = open_file(file_path)

        # Unsupported file type.
        if not mol_dict:
            output_error(msg("invalid_file_format", "sdf", "csv", "json", split=True), pad=1)

    # SDF
    elif ext == "sdf":
        try:
            mol_pdb = Chem.SDMolSupplier(file_path)

            # Check the number of molecules
            molecules = [mol for mol in mol_pdb if mol is not None]
            if len(molecules) == 0:
                raise ValueError("SDF file does not contain valid molecular data")
            if len(molecules) > 1:
                raise ValueError("SDF file contains more than one molecule")

            mol_pdb = molecules[0]
            mol_dict = mol_pdb.GetPropsAsDict()
        except Exception as err:
            output_error(msg("invalid_sdf", err, split=True), pad=1, return_val=False)

    # CSV
    elif ext == "csv":
        import csv

        def _get_delimiter(file_path, bytes=4096):
            sniffer = csv.Sniffer()
            data = open(file_path, "r").read(bytes)
            delimiter = sniffer.sniff(data).delimiter
            return delimiter

        delimiter = _get_delimiter(file_path)

        try:
            df = pd.read_csv(file_path, delimiter=delimiter)

            # Check the number of molecules
            molecules = [mol for mol in df.iterrows() if mol is not None]
            if len(molecules) == 0:
                raise ValueError("CSV file does not contain valid molecular data")
            if len(molecules) > 1:
                raise ValueError("CSV file contains more than one molecule")

            mol_dict = df.to_dict(orient="records")[0]
        except Exception as err:
            output_error(msg("invalid_csv", err, split=True), pad=1, return_val=False)

    # PDB - current molviewer is not equipped to properly display most proteins.
    # Would require different visualization but not a priority since we're focusing
    # on small molecules.
    # elif ext == "pdb":
    #     try:
    #         mol_pdb = Chem.MolFromPDBFile(file_path)

    #         # May need to check for the numbe of molecules.
    #         # This can be done with biopython
    #         mol_dict = {
    #             "InChI": Chem.MolToInchi(mol_pdb),
    #             "SMILES": Chem.MolToSmiles(mol_pdb),
    #         }
    #     except Exception as err:
    #         output_error(msg("invalid_sdf", err, split=True), pad=1)

    return mol_dict


# Create molecule object.
def __create_mol_obj(mol_dict):
    if "InChI" in mol_dict:
        # Prepend 'InChI=' when missing.
        if not mol_dict["InChI"].startswith("InChI="):
            mol_dict["InChI"] = "InChI=" + mol_dict["InChI"]

        return Chem.MolFromInchi(mol_dict["InChI"])
    elif "SMILES" in mol_dict:
        return Chem.MolFromSmiles(mol_dict["SMILES"])


# Fill in missing identifiers: InChI, InChIKey, SMILES.
def _complete_identifiers(mol_dict, mol_obj):
    mol_dict_lowercase = {k.lower(): v for k, v in mol_dict.items()}

    # Add InChI
    if "InChI" not in mol_dict:
        if "inchi" in mol_dict_lowercase:
            mol_dict["InChI"] = mol_dict_lowercase["inchi"]
        else:
            mol_dict["InChI"] = Chem.MolToInchi(mol_obj)

    # Prepend 'InChI=' when missing.
    if not mol_dict["InChI"].startswith("InChI="):
        mol_dict["InChI"] = "InChI=" + mol_dict["InChI"]

    # Add InChIKey
    if "InChIKey" not in mol_dict:
        if "inchikey" in mol_dict_lowercase:
            mol_dict["InChIKey"] = mol_dict_lowercase["inchikey"]
        else:
            mol_dict["InChIKey"] = Chem.inchi.InchiToInchiKey(mol_dict["InChI"])

    # Add SMILES
    if not "SMILES" in mol_dict:
        if "smiles" in mol_dict_lowercase:
            mol_dict["SMILES"] = mol_dict_lowercase["smiles"]
        else:
            mol_dict["SMILES"] = Chem.rdmolfiles.MolToSmiles(mol_obj)


# Create svg code.
def _renderSVG(mol_obj):
    mol_drawer = Chem.Draw.MolDraw2DSVG(300, 300)
    mol_drawer.DrawMolecule(mol_obj)
    mol_drawer.FinishDrawing()
    return mol_drawer.GetDrawingText()


# Create sdf code.
def _mol2sdf(mol_obj):
    # Generate 3D coordinates for the molecule (optional but usually desirable for SDF)
    # Chem.AllChem.EmbedMolecule(mol_obj, Chem.AllChem.ETKDG())

    # Convert molecule object to SDF format
    mol_sdf = Chem.MolToMolBlock(mol_obj)
    return mol_sdf


# # Calculate additional properties.
# def _enrich(mol_dict, mol_obj):
#     props = {
#         #     "MolecularWeight": Chem.rdMolDescriptors.CalcMolWt(mol_obj),
#         "ExactMolecularWeight": Chem.rdMolDescriptors.CalcExactMolWt(mol_obj),
#         # "BCUT2D": Chem.rdMolDescriptors.BCUT2D(mol_obj),
#         # "EEMcharges": Chem.rdMolDescriptors.CalcEEMcharges(mol_obj),
#         # "Eccentricity": Chem.rdMolDescriptors.CalcEccentricity(mol_obj),
#         # "Asphericity": Chem.rdMolDescriptors.CalcAsphericity(mol_obj),
#         # "MolFormula": Chem.rdMolDescriptors.CalcMolFormula(mol_obj),
#         # "NumAliphaticCarbocycles": Chem.rdMolDescriptors.CalcNumAliphaticCarbocycles(mol_obj),
#         # "NumAliphaticHeterocycles": Chem.rdMolDescriptors.CalcNumAliphaticHeterocycles(mol_obj),
#         # "NumAliphaticRings": Chem.rdMolDescriptors.CalcNumAliphaticRings(mol_obj),
#         "NumAmideBonds": Chem.rdMolDescriptors.CalcNumAmideBonds(mol_obj),
#         # "NumAromaticCarbocycles": Chem.rdMolDescriptors.CalcNumAromaticCarbocycles(mol_obj),
#         # "NumAromaticHeterocycles": Chem.rdMolDescriptors.CalcNumAromaticHeterocycles(mol_obj),
#         # "NumAromaticRings": Chem.rdMolDescriptors.CalcNumAromaticRings(mol_obj),
#         "NumAtoms": Chem.rdMolDescriptors.CalcNumAtoms(mol_obj),
#         "NumAtomStereoCenters": Chem.rdMolDescriptors.CalcNumAtomStereoCenters(mol_obj),
#         "NumUnspecifiedAtomStereoCenters": Chem.rdMolDescriptors.CalcNumUnspecifiedAtomStereoCenters(mol_obj),
#         #     # "HeavyAtomCount": Chem.Descriptors.HeavyAtomCount(mol_obj),
#         #     # "HBA": Chem.rdMolDescriptors.CalcNumHBA(mol_obj),  # Number of Hydrogen Bond Acceptors
#         #     # "HBD": Chem.rdMolDescriptors.CalcNumHBD(mol_obj),  # Number of Hydrogen Bond Donors
#         #     # "RotatableBonds": Chem.Descriptors.NumRotatableBonds(mol_obj),
#         #     # "AromaticRings": Chem.rdMolDescriptors.CalcNumAromaticRings(mol_obj),
#         #     # "SaturatedRings": Chem.rdMolDescriptors.CalcNumSaturatedRings(mol_obj),
#         #     # "AliphaticRings": Chem.rdMolDescriptors.CalcNumAliphaticRings(mol_obj),
#         #     # "TPSA": Chem.rdMolDescriptors.CalcTPSA(mol_obj),  # Topological Polar Surface Area
#         #     # "LogP": Chem.Descriptors.MolLogP(mol_obj),
#         #     # "MolarRefractivity": Chem.Descriptors.MolMR(mol_obj),
#         #     # "NumValenceElectrons": Chem.Descriptors.NumValenceElectrons(mol_obj),
#         #     # "NumHeteroatoms": Chem.Lipinski.NumHeteroatoms(mol_obj),
#         #     # "NumRadicalElectrons": Chem.Descriptors.NumRadicalElectrons(mol_obj),
#         #     # "FormalCharge": Chem.rdkit.Chem.rdmolops.GetFormalCharge(mol_obj),
#     }

#     # Merge mol_boj and props
#     mol_dict.update(props)
