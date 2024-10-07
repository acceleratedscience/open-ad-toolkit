"""
Molecules API
"""

import os
import json
import shutil
from urllib.parse import unquote
from flask import Response, request

from openad.smols.smol_functions import (
    get_smol_from_pubchem,
    create_molset_cache_file,
    assemble_cache_path,
    read_molset_from_cache,
    find_smol,
    mws_add,
    mws_remove,
    get_smol_from_mws,
    get_best_available_identifier,
    get_best_available_smiles,
    merge_smols,
)
from openad.smols.smol_transformers import (
    smol2svg,
    smol2mdl,
    molset2dataframe,
    write_dataframe2sdf,
    write_dataframe2csv,
)

# Macromolecule functions
from openad.mmols.mmol_functions import mmol_from_identifier
from openad.mmols.mmol_transformers import mmol2pdb, mmol2cif, cif2mmol

# Helpers
from openad.helpers.files import open_file
from openad.helpers.output import output_error
from openad.helpers.json_decimal_encoder import DecimalEncoder


class MoleculesApi:
    """
    All the API endpoints related to molecules.
    The API endpoints are called from gui_routes.py.
    """

    def __init__(self, cmd_pointer):
        self.cmd_pointer = cmd_pointer

    # -----------------------------
    # Small molecules
    # -----------------------------

    def get_smol_data(self):
        """
        Get molecule data, plus SDF and SVG.
        Used when requesting a molecule by its identifier.
        """

        data = json.loads(request.data) if request.data else {}
        identifier = data["identifier"] if "identifier" in data else ""

        if not identifier:
            response = Response(None, status=500)
            response.status = "No identifier provided."
            return response

        smol = find_smol(self.cmd_pointer, identifier)

        # Fail
        if not smol:
            response = Response(None, status=500)
            response.status = f"No small molecule found with provided identifier '{identifier}'"
            return response

        # Success
        else:
            return smol, 200

    def get_smol_viz_data(self):
        """
        Get a molecule's SVG and SDF data, used to render 2D and 3D visualizations.
        Used when opening a .smol.json file.
        """

        data = json.loads(request.data) if request.data else {}
        inchi_or_smiles = data["inchi_or_smiles"] if "inchi_or_smiles" in data else None

        if not inchi_or_smiles:
            output_error("get_mol_viz_data() -> Missing inchi_or_smiles")
            return "Missing inchi_or_smiles", 500
        try:
            svg = smol2svg(inchi_or_smiles=inchi_or_smiles)
            mdl = smol2mdl(inchi_or_smiles=inchi_or_smiles)

        # When the identifier is not a smiles or inchi
        except Exception:  # pylint: disable=broad-except
            return {
                "error": f"Failed to generate visualisation data from '{inchi_or_smiles}'. Should be inchi or smiles."
            }, 400

        if not svg and not mdl:
            output_error("get_mol_viz_data() -> Failed to generate SVG/MDL from '{inchi_or_smiles}'.")
            return f"Failed to generate visualisation data from '{inchi_or_smiles}'.", 500

        return {"mdl": mdl, "svg": svg}, 200

    def get_mol_data_from_molset(self):
        """
        Get a molecule from a molset file.
        """
        data = json.loads(request.data) if request.data else {}
        cache_id = data["cacheId"] if "cacheId" in data else ""
        index = data["index"] if "index" in data else 1

        # Workspace path
        cache_path = assemble_cache_path(self.cmd_pointer, "molset", cache_id)

        # Read file from cache.
        molset, err_code = open_file(cache_path, return_err="code")

        # Return error
        if err_code:
            return err_code, 500

        return molset[index - 1], 200

    ##

    def add_mol_to_mws(self):
        """
        Add a molecule to the molecule working set.

        Takes either an identifier or a mol object.
        Identifier is slow because the molecule data has to be loaded from PubChem.
        """

        data = json.loads(request.data) if request.data else {}
        smol = data["mol"] if "mol" in data else ""

        # Get best available identifier.
        _, identifier = get_best_available_identifier(smol)

        # Enrich molecule with PubChem data.
        smol_enriched = find_smol(self.cmd_pointer, identifier, basic=True)
        if smol_enriched:
            smol = merge_smols(smol, smol_enriched)

        # Add it to the working set.
        success = mws_add(self.cmd_pointer, smol, force=True)

        return {"status": success}, 200

    def remove_mol_from_mws(self):
        """
        Remove a molecule from your molecule working set.

        Takes either an identifier or a mol object.
        Identifier is slow because the molecule data has to be loaded from PubChem.
        """

        data = json.loads(request.data) if request.data else {}
        smol = data["mol"] if "mol" in data else ""

        # Remove it from the working set.
        success = mws_remove(self.cmd_pointer, smol, force=True)

        return {"status": success}, 200

    def check_mol_in_mws(self):
        """
        Check if a molecule is stored in your molecule working set.
        """

        data = json.loads(request.data) if request.data else {}
        smol = data["mol"] if "mol" in data else ""

        # Get best available identifier.
        _, identifier = get_best_available_identifier(smol)

        # Check if it's in the working set.
        success = bool(get_smol_from_mws(self.cmd_pointer, identifier))

        return {"status": success}, 200

    def enrich_smol(self):
        """
        Enrich a molecule with PubChem data.
        """

        data = json.loads(request.data) if request.data else {}
        smol = data["smol"] if "smol" in data else ""

        # Get best available identifier.
        _, identifier = get_best_available_identifier(smol)

        # Enrich molecule withg PubChem data.
        smol_enriched = get_smol_from_pubchem(identifier)
        if smol_enriched:
            smol = merge_smols(smol, smol_enriched)

        return smol, 200

    ##

    def save_smol_as_json(self):
        """
        Save new .smol.json file to a specified destination path.
        """
        return self._save_mol("mol_json")

    def save_smol_as_sdf(self):
        """
        Save new .sdf file to a specified destination path.
        """
        return self._save_mol("sdf")

    def save_smol_as_csv(self):
        """
        Save new .csv file to a specified destination path.
        """
        return self._save_mol("csv")

    def save_smol_as_mdl(self):
        """
        Save new .mol file to a specified destination path.
        """
        return self._save_mol("mdl")

    def save_smol_as_smiles(self):
        """
        Save new .smi file to a specified destination path.
        """
        return self._save_mol("smiles")

    # -----------------------------
    # Macromolecules
    # -----------------------------

    def get_mmol_data(self):
        """
        Get macromolecule data.
        Used when requesting a macromolecule by its identifier.
        """

        data = json.loads(request.data) if request.data else {}
        identifier = data["identifier"] if "identifier" in data else ""

        if not identifier:
            response = Response(None, status=500)
            response.status = "No identifier provided."
            return response

        success, cif_data = mmol_from_identifier(identifier)

        # Fail
        if not success:
            response = Response(None, status=500)
            response.status = f"No macromolecule found with provided identifier '{identifier}'"
            return response

        # Success
        else:
            mmol = cif2mmol(cif_data)
            return mmol, 200

    ##

    def save_mmol_as_mmol_json(self):
        """
        Save new .pdb.json file to a specified destination path.
        """
        return self._save_mol("mmol_json")

    def save_mmol_as_pdb(self):
        """
        Save new .pdb file to a specified destination path.
        """
        return self._save_mol("pdb")

    def save_mmol_as_cif(self):
        """
        Save new .cif file to a specified destination path.
        """
        return self._save_mol("cif")

    # -----------------------------
    # Molecules shared
    # -----------------------------

    def _save_mol(self, format_as):
        """
        Save a molecule to a file, in the specified format.
        """

        # Note: the new_file parameter is always true for now, but later on
        # when we let users add comments etc, we'll want to be able to update
        # existing files.
        data = json.loads(request.data) if request.data else {}
        new_file = data["newFile"] if "newFile" in data else False
        force = data["force"] if "force" in data else False
        path = unquote(data["path"]) if "path" in data else None
        smol = data["smol"] if "smol" in data else None
        mmol = data["mmol"] if "mmol" in data else None

        if not path:
            return f"save_new_mol() -> Parameter 'path' missing", 500
        if not smol and not mmol:
            return f"save_new_mol() -> Parameter 'smol' or 'mmol' missing", 500

        # Compile path
        workspace_path = self.cmd_pointer.workspace_path()
        file_path = workspace_path + "/" + path

        # Throw error when detination file (does not) exist(s).
        if path:
            if new_file:
                if os.path.exists(file_path) and not force:
                    print(f"_save_mol() -> File already exists: {file_path}")
                    return f"A file with this name already exists.", 409
            else:
                if not os.path.exists(file_path):
                    print(f"_save_mol() -> Destination not found: {file_path}")
                    return f"The file you're trying to save is not found.", 404

        try:
            # -----------------------------
            # Small molecules
            # -----------------------------

            # Save as .smol.json file.
            if format_as == "mol_json":
                # Write to file
                with open(file_path, "w", encoding="utf-8") as f:
                    json.dump(smol, f, ensure_ascii=False, indent=4, cls=DecimalEncoder)

            # Save as .sdf file.
            elif format_as == "sdf":
                df = molset2dataframe([smol], include_romol=True)
                write_dataframe2sdf(df, file_path)

            # Save as .csv file.
            elif format_as == "csv":
                df = molset2dataframe([smol])
                write_dataframe2csv(df, file_path)

            # Save as .mol file.
            elif format_as == "mdl":
                smol2mdl(smol, path=file_path)

            # Save as .smi file.
            elif format_as == "smiles":
                smiles = get_best_available_smiles(smol)
                with open(file_path, "w", encoding="utf-8") as f:
                    f.write(smiles)

            # -----------------------------
            # Mecromolecules
            # -----------------------------

            # Save as .mmol.json file.
            elif format_as == "mmol_json":
                with open(file_path, "w", encoding="utf-8") as f:
                    json.dump(mmol, f, ensure_ascii=False, indent=4, cls=DecimalEncoder)

            # Save as .cif file.
            elif format_as == "cif":
                mmol2cif(mmol, path=file_path)

            # Save as .pdb file.
            elif format_as == "pdb":
                mmol2pdb(mmol, path=file_path)

        # In case the requested file path does not exist.
        # This could only happen if the user changes the folder structure outside
        # of the GUI then tries to save a file without refreshing the browser.
        except FileNotFoundError as err:
            print(f"_save_mol() -> FileNotFoundError: {err} / path: {file_path}")
            return f"The selected destination does not exist.", 404
        except PermissionError:
            print(f"_save_mol() -> PermissionError: {err} / path: {file_path}")
            return f"Access denied", 405
        except Exception as err:
            print(f"_save_mol() -> Exception: {err} / path: {file_path}")
            return f"The selected destination does not exist.", 400

        return "ok", 200

    # -----------------------------
    # Molsets
    # -----------------------------

    def get_molset(self):
        """
        Get a cached molset, filtered by the query.
        Note: opening molset files is handled by fs_attach_file_data() in workers/file_system.py
        """

        data = json.loads(request.data) if request.data else {}
        cache_id = data["cacheId"] if "cacheId" in data else ""
        query = data["query"] if "query" in data else {}

        # Read molset from cache.
        try:
            molset = read_molset_from_cache(self.cmd_pointer, cache_id)
        except ValueError as err:
            return f"get_molset() -> {err}", 500

        # Formulate response object.
        return create_molset_response(molset, query, cache_id), 200

    def get_molset_mws(self):
        """
        Get the list of molecules currently stored in the molecule working set.
        """

        data = json.loads(request.data) if request.data else {}
        query = data["query"] if "query" in data else {}

        if len(self.cmd_pointer.molecule_list) > 0:
            # Compile molset.
            molset = []
            for i, smol in enumerate(self.cmd_pointer.molecule_list):
                smol["index"] = i + 1
                molset.append(smol)

            # Create cache working copy.
            cache_id = create_molset_cache_file(self.cmd_pointer, molset)

            # Read molset from cache.
            try:
                molset = read_molset_from_cache(self.cmd_pointer, cache_id)
            except ValueError as err:
                return f"get_my_mols() -> {err}", 500

            # Formulate response object.
            return create_molset_response(molset, query, cache_id), 200

        else:
            return "empty", 200

    ##

    def remove_from_molset(self):
        """
        Remove molecules from a molset's cached working copy.
        """

        data = json.loads(request.data) if request.data else {}
        cache_id = data["cacheId"] if "cacheId" in data else ""
        indices = data["indices"] if "indices" in data else []
        query = data["query"] if "query" in data else {}

        if not cache_id:
            return f"remove_from_molset() -> Unrecognized cache_id: {cache_id}", 500

        if len(indices) == 0:
            return "remove_from_molset() -> No indices provided", 500

        # Compile path
        cache_path = assemble_cache_path(self.cmd_pointer, "molset", cache_id)

        # Read file from cache
        molset, err_code = open_file(cache_path, return_err="code", dumb=False)  # %%

        # Return error
        if err_code:
            return err_code, 500

        # Remove molecules
        molset = [mol for mol in molset if mol["index"] not in indices]

        # Write to cache
        with open(cache_path, "w", encoding="utf-8") as f:
            json.dump(molset, f, ensure_ascii=False, indent=4, cls=DecimalEncoder)

        # Create response object.
        return create_molset_response(molset, query, cache_id), 200

    def clear_molset_working_copy(self):
        """
        Clear a molset's cached working copy.
        """
        data = json.loads(request.data) if request.data else {}
        cache_id = data["cacheId"] if "cacheId" in data else ""

        if not cache_id:
            return f"clear_molset_working_copy() -> Unrecognized cache_id: {cache_id}", 500

        cache_path = assemble_cache_path(self.cmd_pointer, "molset", cache_id)
        if os.path.exists(cache_path):
            os.remove(cache_path)

        return "ok", 200

    ##

    def update_molset(self):
        """
        Save changes to a molset file.
        """
        return self._save_molset(new_file=False)

    def update_molset_mws(self):
        """
        Save changes to the molecule working set.
        """
        return self._save_molset(format_as="my-mols")

    ##

    def save_molset_as_json(self):
        """
        Save molset as molset.json file.
        """
        data = json.loads(request.data) if request.data else {}
        new_file = data["newFile"] if "newFile" in data else ""
        return self._save_molset(new_file=new_file)

    def save_molset_as_sdf(self):
        """
        Save molset as .sdf file.
        """
        data = json.loads(request.data) if request.data else {}
        remove_invalid_mols = data["removeInvalidMols"] if "removeInvalidMols" in data else ""
        new_file = data["newFile"] if "newFile" in data else ""
        return self._save_molset(new_file=new_file, format_as="sdf", remove_invalid_mols=remove_invalid_mols)

    def save_molset_as_csv(self):
        """
        Save molset as .csv file.
        """
        data = json.loads(request.data) if request.data else {}
        new_file = data["newFile"] if "newFile" in data else ""
        return self._save_molset(new_file=new_file, format_as="csv")

    def save_molset_as_smiles(self):
        """
        Save molset as .smi file.
        """
        data = json.loads(request.data) if request.data else {}
        new_file = data["newFile"] if "newFile" in data else ""
        return self._save_molset(new_file=new_file, format_as="smiles")

    def _save_molset(self, new_file=False, format_as="molset_json", remove_invalid_mols=False):
        """
        Save a molset to a file, in the specified format.

        Parameters
        ----------
        new_file: bool
            Whether we're creating a new file or overwriting an existing one.
        format_as: 'molset_json' | 'sdf' | 'csv' | 'smiles'.
            The format to save the molset as.
        remove_invalid_mols: bool
            Whether to remove invalid molecules from the molset before saving.

        """
        data = json.loads(request.data) if request.data else {}
        cache_id = data["cacheId"] if "cacheId" in data else ""
        path = unquote(data["path"]) if "path" in data else ""

        if not cache_id:
            return f"_save_molset() -> Unrecognized cache_id: {cache_id}", 500

        # Compile path
        workspace_path = self.cmd_pointer.workspace_path()
        file_path = workspace_path + "/" + path
        cache_path = assemble_cache_path(self.cmd_pointer, "molset", cache_id)

        # Throw error when detination file (does not) exist(s).
        if path:
            if new_file:
                if os.path.exists(file_path):
                    return f"_save_molset() -> File already exists: {file_path}", 403
            else:
                if not os.path.exists(file_path):
                    return f"_save_molset() -> File not found: {file_path}", 404

        if not os.path.exists(cache_path):
            return f"_save_molset() -> Cached working copy not found: {cache_path}", 500

        try:
            # For .molset.json files, we simply copy the cache file to the destination.
            if format_as == "molset_json":
                shutil.copy(cache_path, file_path)

            # For all other formats, we need to read the
            # molset data into memory so we can transform it.
            else:
                molset, err = open_file(cache_path, return_err=True)
                if err:
                    return err, 500

            # Save as SDF file.
            if format_as == "sdf":
                try:
                    df = molset2dataframe(molset, remove_invalid_mols, include_romol=True)
                    write_dataframe2sdf(df, file_path)
                except ValueError as err:
                    return {
                        "error": str(err),
                        "invalidMols": err.args[1],
                    }, 422

            # Save as CSV file.
            elif format_as == "csv":
                try:
                    df = molset2dataframe(molset, remove_invalid_mols)
                    write_dataframe2csv(df, file_path)
                except ValueError as err:
                    return {
                        "error": str(err),
                        "invalidMols": err.args[1],
                    }, 422

            # Save as SMILES file.
            elif format_as == "smiles":
                smiles_list = []
                missing_smiles = []
                for mol in molset:
                    smiles = get_best_available_smiles(mol)
                    if smiles:
                        smiles_list.append(smiles)
                    else:
                        missing_smiles.append(mol["index"])

                # Return error if there are missing SMILES.
                if missing_smiles:
                    return {
                        "error": "Some molecules are missing SMILES.",
                        "invalidMols": missing_smiles,
                    }, 422

                # Write to file.
                try:
                    with open(file_path, "w", encoding="utf-8") as f:
                        f.write("\n".join(smiles_list))
                except Exception as err:  # pylint: disable=broad-except
                    return f"Error writing SMILES file: {err}", 500

            elif format_as == "my-mols":
                # Read file from cache.
                molset, err_code = open_file(cache_path, return_err="code")
                if err_code:
                    return err_code, 500

                # Compile molset.
                molecule_list = []
                for mol in molset:
                    try:
                        molecule_list.append(mol)
                    except Exception as fails:  # pylint: disable=broad-except
                        print(f"Error converting molecule format: {fails}", mol)

                self.cmd_pointer.molecule_list = molecule_list

        # In case the requested file path does not exist.
        # This could only happen if the user changes the folder structure outside
        # of the GUI then tries to save a file without refreshing the browser.
        except FileNotFoundError as err:
            return f"_save_molset() -> FileNotFoundError: {err}", 404

        return "ok", 200

    def replace_mol_in_molset(self):
        """
        Replace a molecule in a molset file.

        This first updates the cache working copy, then saves the
        changes to the actual molset file, or to the my-mols list.
        """

        data = json.loads(request.data) if request.data else {}
        cache_id = data["cacheId"] if "cacheId" in data else ""
        mol = data["mol"] if "mol" in data else ""
        context = data["context"] if "context" in data else None
        # TO DO: Sync formatAs in API with context in the GUI
        format_as = "molset_json" if context == "json" else context

        if not cache_id:
            return f"save_mol_in_molset() -> Unrecognized cache_id: {cache_id}", 500

        # Compile path
        cache_path = assemble_cache_path(self.cmd_pointer, "molset", cache_id)

        # Read file from cache.
        molset, err_code = open_file(cache_path, return_err="code")
        if err_code:
            return err_code, 500

        # Replace molecule in molset working copy.
        index = mol.get("index")
        molset[index - 1] = mol

        # Write to cache.
        with open(cache_path, "w", encoding="utf-8") as f:
            json.dump(molset, f, ensure_ascii=False, indent=4, cls=DecimalEncoder)

        # Now the working copy is updated, we also update the molset.
        return self._save_molset(new_file=False, format_as=format_as)


def create_molset_response(molset, query=None, cache_id=None):
    """
    Return a filtered and paginated subset of a molset, wrapped into
    a response object that is ready to be consumed by the frontend.

    Parameters
    ----------
    molset: list
        The entire list of molecules from a file.
    cache_id: str
        The unique identifier of the working copy of the molset.
    query: dict
        The query parameters from the frontend.
        Eg. { queryStr: "C1=CC=CC=C1", smarts: 1, page: 1, sort: "name" }
    """

    # Parse the query
    query = query if query else {}
    search_str = unquote(query["search"]) if "search" in query else ""
    smarts_mode = query["smarts"] if "smarts" in query else False
    page = int(query["page"]) if "page" in query else 1
    sort = query["sort"] if "sort" in query else None
    page_size = 48  # Hardcoded for now

    # Filter by query
    if search_str:
        results = []
        for mol in molset:
            found = False

            # Substructure search - match against smiles only.
            if smarts_mode:
                if search_str.lower() in mol["identifiers"]["canonical_smiles"].lower():
                    results.append(mol)

            # Regular search - match against all identifiers and properties.
            else:
                for key in mol["identifiers"]:
                    if search_str.lower() in str(mol["identifiers"][key]).lower():
                        results.append(mol)
                        found = True
                        break
                if not found:
                    for key in mol["properties"]:
                        if search_str.lower() in str(mol["properties"][key]).lower():
                            results.append(mol)
                            break

    # No query
    else:
        results = molset

    # Sort
    try:
        # print("INDEX:", molset[0].get("index"))
        # print("SORT:", sort)
        if sort:
            reverse = sort.startswith("-")
            sort_key = sort.lstrip("-")
            results = sorted(results, key=lambda mol: _sort_mol(mol, sort_key), reverse=reverse)
    except TypeError as err:
        # In the edge case where our dataset mixes string and
        # number values, we want to avoid crashing the app.
        output_error(err)

    # Store all indices - used by 'select all'
    all_indices = [mol.get("index") for mol in molset]

    # Store matching indices before pagination - used by 'select matching'
    matching_indices = [mol.get("index") for mol in results]

    # Paginate
    total_pages = len(results) // page_size + 1
    page = min(page, total_pages)  # Make sure that page number is lowered in case of a too high value.
    total = len(molset)
    result_count = len(results)
    skip = 48 * (page - 1)
    results = results[skip : skip + 48]

    return {
        "cacheId": cache_id,  # Used to identify our working copy in next requests
        "mols": results,  # One page of molecules
        "allIndices": all_indices,  # Ids of all molecules
        "matchingIndices": matching_indices,  # Ids of all matching molecules
        "total": total,
        "resultCount": result_count,
        # Query parameters:
        "searchStr": search_str,
        "searchMode": "smarts" if smarts_mode else "text",
        "sort": sort,
        "page": page,
        "pageSize": page_size,
    }


def _sort_mol(mol, sort_key):
    """
    Sorter function for a molset.

    Parameters
    ----------
    mol: dict
        A molecule object.
    sort_key: str
        The key of the category whose value we'll sort by.
        Eg. 'name' (identifier) or 'molecular_weight' (property).
    """
    if sort_key == "index":
        value = mol.get(sort_key)
    elif sort_key == "name":
        value = (mol.get("identifiers") or {}).get(sort_key)
    else:
        value = (mol.get("properties") or {}).get(sort_key)

    value = __prep_sort_value(value)

    # Returning a tuple will sort by the first value, then the
    # second, etc. This lets us group all none values on top.
    return (value is None, value)


def __prep_sort_value(value):
    """
    Prepare a value for sorting by:
    - Converting strings to lowercase
    - Converting number strings to floats

    Parameters
    ----------
    value: str, int, float
        The value to prepare.
    """

    # Convert number strings to floats
    try:
        return float(value)
    except ValueError:
        pass

    # Convert strings to lowercase
    if isinstance(value, str):
        return value.lower()

    return value
