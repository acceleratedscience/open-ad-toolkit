# MAJOR-RELEASE-TODO:
# Certain commands are Jupyter-only (or should be separated as jupyter only) and should be hidden from CLI help.
# Examples:
# - there should be a `fetch molecule` command, or `@molecule_x` that returns molecule data
# - The `export molecule|mol` / `export molecules|mols` comands returns molecule data but this should be replaced with the command above.
#   Other command like `display data` do the same, adding to the inconsistency and confusion.


from openad.core.help import help_dict_create
from openad.app.global_var_lib import GLOBAL_SETTINGS
from openad.smols.smol_functions import SMOL_PROPERTIES, PCY_IDFR
from openad.helpers.general import is_notebook_mode
from openad.helpers.pretty_data import list_columns

from pyparsing import (
    alphanums,
    alphas,
    CaselessKeyword,
    CharsNotIn,
    Combine,
    delimitedList,
    Forward,
    Group,
    Keyword,
    MatchFirst,
    nums,
    OneOrMore,
    Optional,
    ParserElement,
    QuotedString,
    Suppress,
    Word,
    ZeroOrMore,
    # Literal,
    # replaceWith,
    # Combine,
    # pyparsing_test,
    # ParseException,
)

(
    add,
    append,
    a_s,
    basic,
    context,
    create,
    create,
    data,
    description,
    d_isplay,
    e_xec,
    experiment,
    export,
    file,
    force,
    formula,
    f_rom,
    get,
    GPU,
    history,
    inchi,
    inchikey,
    jobs,
    l_ast,
    l_ist,
    load,
    merge,
    name,
    only,
    optimize,
    protein,
    pubchem,
    remove,
    rename,
    result,
    results,
    run,
    runs,
    save,
    s_et,
    show,
    smiles,
    sources,
    toolkit,
    toolkits,
    unset,
    upsert,
    using,
    w_ith,
    workspace,
    workspaces,
) = map(
    CaselessKeyword,
    "add append as basic context create create data description display exec experiment export file force formula from get \
    gpu history inchi inchikey jobs last list load merge name only optimize protein pubchem remove rename result results \
    run runs save set show smiles sources toolkit toolkits unset upsert using with workspace workspaces".split(),
)
_mol = ["molecule", "mol"]
_mols = ["molecules", "mols"]
_molset = ["molecule-set", "molset"]
_molsets = ["molecule-sets", "molsets"]
_mol_properties = ["synonyms"]
_mol_properties.extend(SMOL_PROPERTIES)
clear = CaselessKeyword("clear")
cache = CaselessKeyword("cache")
analysis = CaselessKeyword("analysis")
enrich = CaselessKeyword("enrich")
mol_properties = MatchFirst(map(CaselessKeyword, _mol_properties))
molecules = MatchFirst(map(CaselessKeyword, _mols))
molecule = MatchFirst(map(CaselessKeyword, _mol))
molecule_set = MatchFirst(map(CaselessKeyword, _molset))
molecule_sets = MatchFirst(map(CaselessKeyword, _molsets))
molecule_identifier = Word(
    alphas, alphanums + "_" + "[" + "]" + "(" + ")" + "=" + "," + "-" + "+" + "/" + "#" + "@" + "." + "*" + ";"
) | Word(nums)
protein_identifier = Word(alphas + "-", "*")
desc = QuotedString("'", escQuote="\\")

MOL_SHORTHAND = "You can use the 'mol' shorthand instead of 'molecule'."
MOLS_SHORTHAND = "You can use the 'mols' shorthand instead of 'molecules'."
MOLSET_SHORTHAND = "You can use the 'molset' shorthand instead of 'molecule-set'."
MOL_LOOKUP_PRIORITY = (
    "If the requested molecule exists in your current working set or in memory, that version will be prioritized."
)
SUPPORTED_IDENTIFIERS = """Supported molecule identifiers:
- <cmd>name</cmd> / <cmd>synonym</cmd>
- <cmd>SMILES</cmd>
- <cmd>InChI</cmd>
- <cmd>InChIKey</cmd>
- <cmd>PubChem CID</cmd>"""
SUPPORTED_IDENTIFIERS_BASIC = """Supported molecule identifiers:
- <cmd>name</cmd> / <cmd>synonym</cmd>
- <cmd>SMILES</cmd> <soft>- supports [ basic ]</soft>
- <cmd>InChI</cmd> <soft>- supports [ basic ]</soft>
- <cmd>InChIKey</cmd>
- <cmd>PubChem CID</cmd>"""
SUPPORTED_FILE_FORMATS = """Supported file formats:
- molset (.molset.json)
- SDF (.sdf)
- CSV (.csv)
- SMILES (.smi)"""
EXAMPLE_INPUT_FILES = "To see some example input files, you can export the molecules from your working set using the <cmd>export molecules|mols as ...</cmd> command."


def smol_grammar_add(statements, grammar_help):
    """
    Grammar for managing small molecules.
    """

    #
    #
    # Small Molecules
    #
    #

    # ---
    # Display a molecule's details
    statements.append(
        Forward(d_isplay + molecule + (molecule_identifier | desc)("molecule_identifier"))("display_molecule")
    )
    grammar_help.append(
        help_dict_create(
            name="display molecule",
            category="Small Molecules",
            command="display molecule|mol <name> | <smiles> | <inchi> | <inchikey> | <cid>",
            description=f"""Display a molecule's details.

A molecule's details include its identifiers, synonyms, properties and any analysis results it has been enriched with.

{SUPPORTED_IDENTIFIERS}

Notes:
- {MOL_SHORTHAND}
- {MOL_LOOKUP_PRIORITY}
            
Examples:
- Display a molecule by name:
  <cmd>display molecule Aspirin</cmd>

- Display a molecule by SMILES string:
  <cmd>display molecule CC(=O)OC1=CC=CC=C1C(=O)O</cmd>

- Display a molecule by InChI string:
  <cmd>display mol InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)</cmd>

- Display a molecule by InChIKey:
  <cmd>display mol BSYNRYMUTXBXSQ-UHFFFAOYSA-N</cmd>

- Display a molecule by its PubChem CID:
  <cmd>display mol 2244</cmd>
""",
        )
    )

    # ---
    # Show individual molecule in browser
    statements.append(
        Forward(show("show") + molecule + (molecule_identifier | desc)("molecule_identifier"))("show_mol")
    )  # From mol json file
    grammar_help.append(
        help_dict_create(
            name="show molecule",
            category="Small Molecules",
            command="show molecule|mol <name> | <smiles> | <inchi> | <inchikey> | <cid>",
            description=f"""Launch the molecule viewer { 'in your browser ' if is_notebook_mode() else '' }to visualize your molecule and inspect its properties.

{SUPPORTED_IDENTIFIERS}

Notes:
- {MOL_SHORTHAND}
- {MOL_LOOKUP_PRIORITY}

Examples:
- Inspect a molecule by name:
  <cmd>show mol aspirin</cmd>

- Inspect a molecule by SMILES string:
  <cmd>show mol CC(=O)OC1=CC=CC=C1C(=O)O</cmd>

- Inspect a molecule by InChI string:
  <cmd>show mol InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)</cmd>

- Inspect a molecule by PubChem CID:
  <cmd>show mol 2244</cmd>
""",
        )
    )

    # ---
    # Show molset in browser from file / dataframe
    statements.append(Forward(show("show") + molecule_set + desc("molset_file"))("show_molset"))
    statements.append(
        Forward(
            show("show")
            + molecule_set
            + using
            + CaselessKeyword("dataframe")
            + Word(alphas, alphanums + "_")("in_dataframe")
        )("show_molset_df")
    )
    grammar_help.append(
        help_dict_create(
            name="show molset",
            category="Small Molecules",
            command="show molset|molecule-set '<molset_or_sdf_or_smi_path>' | using dataframe <dataframe>",
            description=f"""Launch the molset viewer { 'in your browser ' if is_notebook_mode() else '' }to visualize the contents of a molecule set.

Examples:
- <cmd>show molset 'neurotransmitters.mol.json'</cmd>
- <cmd>show molset 'neurotransmitters.sdf'</cmd>
- <cmd>show molset 'neurotransmitters.smi'</cmd>
- <cmd>show molset my_dataframe</cmd>
""",
        )
    )

    # ---
    # Get molecule property
    _mol_properties = ["synonyms"]

    _mol_properties.extend(SMOL_PROPERTIES)
    _mol_properties.extend(PCY_IDFR.keys())
    mol_properties = MatchFirst(map(Keyword, _mol_properties))

    statements.append(
        Forward("@" + (molecule_identifier | desc)("molecule_identifier") + ">>" + mol_properties("property"))(
            "get_smol_prop"
        )
    )
    # statements.append(
    #    Forward("@" + (molecule_identifier | desc)("molecule_identifier") + ">>" + Word(alphas, alphanums + "_"))(
    #        "get_smol_prop_lookup_error"
    #    )
    # )
    grammar_help.append(
        help_dict_create(
            name="@<molecule>",
            category="Small Molecules",
            command="@(<name> | <smiles> | <inchi> | <inchikey> | <cid>)>><molecule_property_name>",
            description=f"""Request a molecule's property.

{SUPPORTED_IDENTIFIERS}

Notes:
- In addition to a molecule's properties or identifiers, you can also request its synonyms.
- {MOL_LOOKUP_PRIORITY}
- In addition to the properties listed below, you can also request any additional properties that are available in your molecule working set.
            
Examples:
- Obtain the molecular weight of a molecule known as Aspirin.
  <cmd>@aspirin>>molecular_weight</cmd>

- Obtain the canonical smiles string for the dopamine neurotransmitter.
  <cmd>@dopamine>>canonical_smiles</cmd>

- Obtain a molecules xlogp value using a SMILES string.
  <cmd>@CC(=O)OC1=CC=CC=C1C(=O)O>>xlogp</cmd>

- Obtain all the synonyms for ibuprofen.
 <cmd>@ibuprofen>>synonyms</cmd>

Available properties that can be queried:
{list_columns(_mol_properties)}
""",
        )
    )

    statements.append(
        Forward(
            export
            + molecule
            + (molecule_identifier | desc)("molecule_identifier")
            + Optional(CaselessKeyword("as") + CaselessKeyword("file"))("as_file")
        )("export_molecule")
    )
    # @later - the way this works is hard to explain and not very intuitive. I would suggest making "as file" mandatory in CLI, so behavior is consistent for both CLI and Notebook.
    grammar_help.append(
        help_dict_create(
            name="export molecule",
            category="Small Molecules",
            command="export molecule|mol <name> | <smiles> | <inchi> | <inchikey> | <cid> [ as file ]",
            description=f"""Either save a molecule to your workspace as a JSON file (CLI) or return a dictionary of the molecule's properties (Jupyter Notebook).

{SUPPORTED_IDENTIFIERS}
            
Notes:
- The requested molecule does not need to be in your current working set.
- The <cmd>as file</cmd> clause is only needed when you wish to save a file to your workspace from within a Jupytyer Notebook.
- {MOL_SHORTHAND}
- {MOL_LOOKUP_PRIORITY}

Examples:
- Export a molecule by name:
  <cmd>export molecule aspirin</cmd>

- Export a molecule by SMILES string, and force it to save as a file (only relevant in a Jupyter Notebook):
  <cmd>export mol aspirin as file</cmd>
""",
        )
    )

    #
    #
    # Molecule Working Set
    #
    #

    # ---
    # Add molecule to working set
    statements.append(
        Forward(
            add
            + molecule
            + (molecule_identifier | desc)("molecule_identifier")
            + Optional(a_s + desc("name"))
            + Optional(basic)("basic")
            + Optional(force)("force")
        )("add_molecule")
    )
    grammar_help.append(
        help_dict_create(
            name="add molecule",
            category="Molecule Working Set",
            command="add molecule|mol <name> | <smiles> | <inchi> | <inchikey> | <cid> [ as '<name>' ] [ basic ] [ force ]",
            description=f"""Add a molecule to your current molecule working set.

{SUPPORTED_IDENTIFIERS_BASIC}

Options:
- <cmd>as <name></cmd>: Provide a custom name for the molecule, which will be used by the software whenever refering to it going forward.
  Note: you can always update a molecule's name later by running <cmd>rename molecule <name></cmd>.
- <cmd>basic</cmd>: Create a minimal molecule without enriching it with PubChem data. This is only relevant when using a SMILES or InChI string as identifier. Because no API calls are made, this is much faster than the default behavior.
- <cmd>force</cmd>: Suppress the confirmation step before adding a molecule, which may be desired in batch operations.

Notes:
- {MOL_SHORTHAND}
- {MOL_LOOKUP_PRIORITY}

Examples:
- Add a molecule by SMILES string:
  <cmd>add molecule CC(=O)OC1=CC=CC=C1C(=O)O</cmd>

- Add a molecule by SMILES string, without enriching it with PubChem data:
  <cmd>add molecule CC(=O)OC1=CC=CC=C1C(=O)O basic</cmd>

- Add a molecule by SMILES string, giving it a custom name:
  <cmd>add molecule CC(=O)OC1=CC=CC=C1C(=O)O as 'mymol' basic</cmd>

- Add a molecule by unquoted InChI string:
  <cmd>add mol InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)</cmd>

- Add a molecule by quoted InChI string:
  <cmd>add mol 'InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)'</cmd>

- Add a molecule by InChIKey:
  <cmd>add mol BSYNRYMUTXBXSQ-UHFFFAOYSA-N</cmd>

- Add a molecule by InChIKey, giving it a custom name:
  <cmd>add mol BSYNRYMUTXBXSQ-UHFFFAOYSA-N as 'mymol'</cmd>

- Add a molecule by name:
  <cmd>add molecule aspirin</cmd>

- Add a molecule by name, supressing the confirmation step:
  <cmd>add molecule aspirin force</cmd>

- Add a molecule by quoted name:
  <cmd>add mol 'Aspirin 325 mg'</cmd>

- Add a molecule by its PubChem CID:
  <cmd>add mol 2244</cmd>
""",
        )
    )

    # ---
    # Remove molecule from working set
    statements.append(
        Forward(remove + molecule + (molecule_identifier | desc)("molecule_identifier") + Optional(force)("force"))(
            "remove_molecule"
        )
    )
    grammar_help.append(
        help_dict_create(
            name="remove molecule",
            category="Molecule Working Set",
            command="remove molecule|mol <name> | <smiles> | <inchi> | <inchikey> | <cid> [ force ]",
            description=f"""Remove a molecule from the current working set based on a given molecule identifier.

{SUPPORTED_IDENTIFIERS}

Options:
- <cmd>force</cmd>: Suppress the confirmation step before removing a molecule, which may be desired in batch operations.

Notes:
- {MOL_SHORTHAND}
            
Examples:
- Remove a molecule by name:
  <cmd>remove molecule Aspirin</cmd>
  
- Remove a molecule by SMILES:
  <cmd>remove molecule CC(=O)OC1=CC=CC=C1C(=O)O</cmd>

- Remove a molecule by InChIKey:
  <cmd>remove mol BSYNRYMUTXBXSQ-UHFFFAOYSA-N</cmd>

- Remove a molecule by InChI
  <cmd>remove mol  InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)</cmd>

- Remove a molecule by CID:
  <cmd>remove mol 2244</cmd>
""",
        )
    )

    # ---
    # List molecules in working set
    statements.append(Forward(l_ist + molecules)("list_molecules"))
    grammar_help.append(
        help_dict_create(
            name="list molecules",
            category="Molecule Working Set",
            command="list molecules|mols",
            description=f"""List all molecules in the current working set.

Notes:
- {MOLS_SHORTHAND}
""",
        )
    )

    # ---
    # Show molecule working set in browser/iFrame (CLI/Jupyter)
    statements.append(Forward(show("show") + molecules)("show_molecules"))
    grammar_help.append(
        help_dict_create(
            name="show molecules",
            category="Molecule Working Set",
            command="show molecules|mols",
            description=f"""Launch the molset viewer { 'in your browser ' if is_notebook_mode() else '' }to visualize your molecule working set.

Notes:
- {MOLS_SHORTHAND}
""",
        )
    )

    # ---
    # Enrich molecule set
    statements.append(Forward(enrich + molecules + w_ith + analysis)("enrich_mws_with_analysis"))
    grammar_help.append(
        help_dict_create(
            name="enrich molecules",
            category="Molecule Working Set",
            command="enrich molecules|mols with analysis",
            description=f"""Enrich the molecules in your current working set with the results of the last performed analysis.

This assumes that your molecule working set contains either the input molecule or any of the result molecules from the analysis.

Notes:
- {MOLS_SHORTHAND}

Currently supported analysis commands:

RXN:
- <cmd>predict reaction</cmd>
- <cmd>predict retrosynthesis</cmd>

DS4SD:
- <cmd>search for patents containing molecule</cmd>
- <cmd>search for similiar molecules</cmd>

Please refer to the DS4SD and RXN toolkits for further assistance on these commands. For example:
<cmd>set context rxn</cmd>
<cmd>predict reaction ?</cmd>
""",
        )
    )

    # ---
    # Clear analysis cache
    statements.append(Forward(clear + analysis + cache)("clear_analysis"))
    grammar_help.append(
        help_dict_create(
            name="clear analysis cache",
            category="Molecule Working Set",
            command="clear analysis cache",
            description="""Clear the analysis results cache for your current workspace.

Please refer to the <cmd>enrich molecules|mols with analysis</cmd> command for more information about analysis results.
""",
        )
    )

    # ---
    # Display working set molecule's sources
    statements.append(
        Forward(d_isplay + sources + (molecule_identifier | desc)("molecule_identifier"))("display_property_sources")
    )
    grammar_help.append(
        help_dict_create(
            name="display sources",
            category="Molecule Working Set",
            command="display sources <name> | <smiles> | <inchi> | <inchikey> | <cid>",
            description=f"""Display the sources of a molecule's properties, attributing how they were calculated or where they were sourced.

{SUPPORTED_IDENTIFIERS}

Notes:
- {MOL_LOOKUP_PRIORITY}
""",
        )
    )

    # ---
    # Rename molecule
    statements.append(
        Forward(
            rename
            + molecule
            + (molecule_identifier | desc)("molecule_identifier")
            + a_s
            + (Word(alphas, alphanums + "_")("new_name"))
        )("rename_molecule")
    )
    grammar_help.append(
        help_dict_create(
            name="rename molecule",
            category="Molecule Working Set",
            command="rename molecule|mol <molecule_identifer_string> as <molecule_name>",
            description=f"""Rename a molecule in the current working set.

Notes:
- {MOL_SHORTHAND}

Example:
- Assuming you've added the molecule <yellow>CC(=O)OC1=CC=CC=C1C(=O)O</yellow> to your molecule working set, you can then rename it as such:
  <cmd>rename molecule CC(=O)OC1=CC=CC=C1C(=O)O as Aspirin</cmd>
""",
        )
    )

    # ---
    # Load molecules from file
    statements.append(
        Forward(
            load
            + molecules
            + f_rom
            + file
            + desc("moles_file")
            + Optional(enrich)("enrich_pubchem")
            + Optional(CaselessKeyword("append"))("append")
        )("load_molecules_file")
    )
    # DEPRECATED: Backward compatibility
    statements.append(
        Forward(
            load
            + molecules
            + using  # <-- changed
            + file
            + desc("moles_file")
            + Optional((merge + w_ith + pubchem))("enrich_pubchem")  # <-- changed
            + Optional(CaselessKeyword("append"))("append")
        )("load_molecules_file-DEPRECATED")
    )
    grammar_help.append(
        help_dict_create(
            name="load molecules",
            category="Molecule Working Set",
            command="load molecules|mols from file '<filename.molset.json|sdf|csv|smi>' [ enrich ] [ append ]",
            description=f"""Load molecules from a file into your molecule working set.

{SUPPORTED_FILE_FORMATS}

Options:
- Append <cmd>enrich</cmd> to enrich the molecule with data from pubchem.
- Append <cmd>append</cmd> to append the molecules to the existing working set instead of overwriting it.

Notes:
- {EXAMPLE_INPUT_FILES}
- {MOLS_SHORTHAND}

Examples:
- Load molecules from a molset JSON file:
  <cmd>load molecules from file 'my_molecules.molset.json'</cmd>

- Load molecules from an SDF file, appending them to the existing working set:
  <cmd>load mols from file 'my_molecules.sdf'` append</cmd>

- Load molecules from an CSV file:
  <cmd>load molecules from file 'my_molecules.csv'</cmd>

- Load molecules from an SMILES file, enriching them with PubChem data:
  <cmd>load mols from file 'my_molecules.smi'` enrich</cmd>
""",
        )
    )

    # ---
    # Load molecules from dataframe
    statements.append(
        Forward(
            load
            + molecules
            + f_rom
            + CaselessKeyword("dataframe")
            + molecule_identifier("in_dataframe")
            + Optional(enrich)("enrich_pubchem")
            + Optional(CaselessKeyword("append"))("append")
        )("load_molecules_dataframe")
    )
    # DEPRECATED: Backward compatibility
    statements.append(
        Forward(
            load
            + molecules
            + using  # <-- changed
            + CaselessKeyword("dataframe")
            + molecule_identifier("in_dataframe")
            + Optional((merge + w_ith + pubchem))("enrich_pubchem")
            + Optional(CaselessKeyword("append"))("append")
        )("load_molecules_dataframe-DEPRECATED")
    )
    grammar_help.append(
        help_dict_create(
            name="load Molecule Working Set",
            category="Molecule Working Set",
            command="load molecules|mols from dataframe <dataframe> [ enrich ] [ append ]",
            description=f"""Load molecules from a dataframe into your molecule working set.

Options:
- Append <cmd>enrich</cmd> to enrich the molecule with data from pubchem.
- Append <cmd>append</cmd> to append the molecules to the existing working set instead of overwriting it.

Notes:
- This command only works when called from a Jupyter Notebook or the API.
- {EXAMPLE_INPUT_FILES}
- {MOLS_SHORTHAND}

Examples:
- Load molecules from a dataframe:
  <cmd>load molecules from dataframe my_dataframe</cmd>

- Load molecules from a dataframe, appending them to the existing working set:
  <cmd>load mols from dataframe my_dataframe append</cmd>

- Load molecules from a dataframe, enriching them with PubChem data:
  <cmd>load mols from dataframe my_dataframe enrich</cmd>
""",
        )
    )

    # ---
    # Merge molecules from dataframe
    statements.append(
        Forward(
            CaselessKeyword("merge")
            + molecules
            + data
            + f_rom
            + CaselessKeyword("dataframe")
            + molecule_identifier("in_dataframe")
            + Optional(enrich)("enrich_pubchem")
        )("merge_molecules_data_dataframe")
    )
    # DEPRECATED: Backward compatibility
    statements.append(
        Forward(
            CaselessKeyword("merge")
            + molecules
            + data
            + using  # <-- changed
            + CaselessKeyword("dataframe")
            + molecule_identifier("in_dataframe")
            + Optional((merge + w_ith + pubchem))("enrich_pubchem")  # <-- changed
        )("merge_molecules_data_dataframe-DEPRECATED")
    )
    grammar_help.append(
        help_dict_create(
            name="merge molecules data",
            category="Molecule Working Set",
            command="merge molecules|mols data from dataframe <dataframe> [ enrich ]",
            description=f"""Merge molecule data from a dataframe into the molecules in your working set.

Options:
- Append <cmd>enrich</cmd> to enrich the molecule with data from pubchem.

The dataframe columns should be named as follows:
- <cmd>subject</cmd> or <cmd>smiles</cmd>: molecules similes string
- <cmd>property</cmd>: the name of the property to be merged
- <cmd>result</cmd>: the value of the property to be nmerged

Sample input file:

<reverse> subject                                                               property                        result </reverse>
<reverse> --------------------------------------------------------------------  -------------------------  ----------- </reverse>
<reverse> O=C(N)C(F)(OC(F)(F)C(F)(F)C(F)(F)F)C(F)(F)F                           molecular_weight               329.065 </reverse>
<reverse> O=C(O)C(F)(NC(F)(F)C(F)(F)C(F)(F)F)C(F)(F)F                           molecular_weight               329.065 </reverse>
<reverse> O=C(O)C(F)(OC(F)(F)C(F)(F)C)CF                                        molecular_weight               240.099 </reverse>
<reverse> O=C(O)C(F)(OC(O)(F)C(F)(F)C(F)(F)F)C(F)(F)F                           molecular_weight               328.058 </reverse>
<reverse> O=C(O)C(F)(OC(Cl)(F)C(F)(F)C(F)(F)F)C(F)(F)F                          molecular_weight               346.504 </reverse>
<reverse> O=C(O)C(F)(OC(F)(F)C(F)(O)C(F)(F)F)C(F)(F)F                           molecular_weight               328.058 </reverse>
<reverse> O=C(O)C(F)(OC(F)(O)C(F)(F)C(F)(F)F)C(F)(F)F                           molecular_weight               328.058 </reverse>
<reverse> O=C(O)C(F)(OC(F)(F)C(F)(Br)C(F)(F)F)C(F)(F)F                          molecular_weight               390.955 </reverse>
<reverse> O=C(O)C(F)OC(O)(F)C(F)(F)C(F)(F)F                                     molecular_weight               260.061 </reverse>

Notes:
- {MOL_SHORTHAND}

Examples:
- Merge molecule data from a dataframe called <cmd>new_props</cmd>:
  <cmd>merge molecules data from dataframe new_props</cmd>

- Merge molecule data from a dataframe called <cmd>new_props</cmd>, while enriching the molecules with PubChem data:
  <cmd>merge molecules data from dataframe new_props enrich</cmd>
""",
        )
    )

    # ---
    # Export all molecules ferom working set
    statements.append(Forward((export + molecules + Optional(a_s + desc("file_name")))("export_mws")))
    grammar_help.append(
        help_dict_create(
            name="export molecules",
            category="Molecule Working Set",
            command="export molecules|mols [ as '<filename.molset.json|sdf|csv|smi>' ]",
            description=f"""Export your molecule working set as a file (CLI) or return it as a dataframe (Jupyter/API).

{SUPPORTED_FILE_FORMATS}

Notes:
- When exporting as a file, the filename's extension will define what format the molecule are exported as.
- If no filename or extension is provided, the molecules will be saved as CSV file.
- When run inside a Jupyter Notebook or from the API, the <cmd>as <filename></cmd> clause will be ignored and a dataframe will be returned.
- {MOLS_SHORTHAND}
""",
        )
    )

    # ---
    # Clear molecule working set
    statements.append(Forward(clear + molecules + Optional(force)("force"))("clear_molecules"))
    grammar_help.append(
        help_dict_create(
            name="clear Molecules",
            category="Molecule Working Set",
            command="clear molecules|mols [ force ]",
            description=f"""Clear the molecule working set.

Options:
- <cmd>force</cmd>: Suppress the confirmation step before clearing the working set, which may be desired in batch operations.

Notes:
- {MOLS_SHORTHAND}
""",
        )
    )

    #
    #
    # DEPRECATED - Remove at next major release / MAJOR-RELEASE-TODO
    #
    #

    # Export molecules as .molecule files
    # $ save molset as foobar
    statements.append(
        Forward(save + molecule_set + a_s + Word(alphas, alphanums + "_")("molset_name"))("save_molecules_DEPRECATED")
    )

    # Import molecules from .molecule files
    # $ load molset foobar
    statements.append(
        Forward(load + molecule_set + Word(alphas, alphanums + "_")("molset_name"))("load_molecules_DEPRECATED")
    )

    # Merge molecules from .molecule files
    # $ merge molset foobar
    # $ merge molset foobar merge only
    # fmt: off
    statements.append(Forward(merge + molecule_set + Word(alphas, alphanums + "_")("molset_name") + Optional(merge + only)("merge_only") + Optional(append + only)("append_only"))("merge_molecules_DEPRECATED"))
    # fmt: on

    # List old molecule sets (folders of .molecule files) in your workspace
    # $ list molsets
    statements.append(Forward(l_ist + molecule_sets)("list_molecule_sets_DEPRECATED"))
