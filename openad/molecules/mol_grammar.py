from openad.core.help import help_dict_create
from openad.molecules.mol_functions import MOL_PROPERTIES as m_props
from openad.helpers.general import is_notebook_mode

m_props.sort()
from pyparsing import (
    Word,
    delimitedList,
    alphas,
    alphanums,
    OneOrMore,
    ZeroOrMore,
    CharsNotIn,
    Forward,
    CaselessKeyword,
    MatchFirst,
    Keyword,
    QuotedString,
    ParserElement,
    Suppress,
    Optional,
    Group,
    Combine,
    nums,
    # Literal,
    # replaceWith,
    # Combine,
    # pyparsing_test,
    # ParseException,
)

(
    get,
    lister,
    description,
    using,
    create,
    s_et,
    unset,
    workspace,
    workspaces,
    context,
    jobs,
    e_xec,
    a_s,
    optimize,
    w_ith,
    toolkits,
    toolkit,
    GPU,
    experiment,
    add,
    run,
    save,
    runs,
    show,
    file,
    d_isplay,
    history,
    data,
    remove,
    result,
    f_rom,
    inchi,
    inchikey,
    smiles,
    formula,
    name,
    l_ast,
    load,
    results,
    export,
    create,
    rename,
    merge,
    pubchem,
    sources,
) = map(
    CaselessKeyword,
    "get list description using create set unset workspace workspaces context jobs exec\
          as optimize with toolkits toolkit gpu experiment add run save runs show \
              file display history data remove result from inchi inchikey smiles formula name last load results export create rename merge pubchem sources".split(),
)
mol = ["molecule", "mol"]
mols = ["molecules", "mols"]
molset = ["molecule-set", "molset"]
molsets = ["molecule-sets", "molsets"]
clear = CaselessKeyword("clear")
cache = CaselessKeyword("cache")
analysis = CaselessKeyword("analysis")
enrich = CaselessKeyword("enrich")
mol_properties = ["synonyms"]
mol_properties.extend(m_props)
mol_properties = MatchFirst(map(CaselessKeyword, mol_properties))
molecules = MatchFirst(map(CaselessKeyword, mols))
molecule = MatchFirst(map(CaselessKeyword, mol))
molecule_set = MatchFirst(map(CaselessKeyword, molset))
molecule_sets = MatchFirst(map(CaselessKeyword, molsets))
molecule_identifier = Word(
    alphas, alphanums + "_" + "[" + "]" + "(" + ")" + "=" + "," + "-" + "+" + "/" + "#" + "@" + "." + "*" + ";"
) | Word(nums)


desc = QuotedString("'", escQuote="\\")

MOL_SHORTHAND = "You can use the 'mol' shorthand instead of 'molecule'."
MOLSET_SHORTHAND = "You can use the 'molset' shorthand instead of 'molecule-set'."
SPECIFY_MOL = "You can specify any molecule by SMILES or InChI, and PubChem classified molecules also by name, InChIKey or their PubChem CID."
WORKING_SET_PRIORITY = "If the requested molecule exists in your current working set, that version will be used."
USING_NAME = "If you use the name of a molecule, the tool will do a caseless search of the names and synonyms first in current working set, then on pubchem."


def mol_grammar_add(statements, grammar_help):
    """defines the grammar available for managing molecules"""

    # ---
    # Add molecule
    statements.append(Forward(add + molecule + molecule_identifier("molecule_identifier"))("add_molecule"))
    grammar_help.append(
        help_dict_create(
            name="add molecule",
            category="Molecules",
            command="add molecule <name> | <smiles> | <inchi> | <inchikey> | <cid>",
            description=f"""
Add a molecule to the current working set of molecules.

{MOL_SHORTHAND}

{SPECIFY_MOL}

{USING_NAME}

When adding a molecule by name, this name will become the molecule's identifying string. You can set or override an identifying string for any molecule with the <cmd>rename molecule</cmd> command.
            
Examples:
- Add a molecule by name:
<cmd>add molecule aspirin</cmd>

- Add a molecule by SMILES:
<cmd>add molecule CC(=O)OC1=CC=CC=C1C(=O)O</cmd>

- Add a molecule by CID:
<cmd>add mol 2244</cmd>

- Add a molecule by InChI:
<cmd>add mol InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)</cmd>

- Add a molecule by InChIKey:
<cmd>add mol BSYNRYMUTXBXSQ-UHFFFAOYSA-N</cmd>
""",
        )
    )

    # ---
    # Display molecule
    statements.append(Forward(d_isplay + molecule + (molecule_identifier)("molecule_identifier"))("display_molecule"))
    grammar_help.append(
        help_dict_create(
            name="display molecule",
            category="Molecules",
            command="display molecule <name> | <smiles> | <inchi> | <inchikey> | <cid>",
            description=f"""
Display a molecule's properties.

{MOL_SHORTHAND}

{WORKING_SET_PRIORITY}

{SPECIFY_MOL}

{USING_NAME}
            
Examples:
- Display a molecule by name:
<cmd>display molecule Aspirin</cmd>

- Display a molecule by SMILES:
<cmd>display molecule CC(=O)OC1=CC=CC=C1C(=O)O</cmd>

- Display a molecule by InChI:
<cmd>display mol InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)</cmd>

- Display a molecule by InChIKey string:
<cmd>display mol BSYNRYMUTXBXSQ-UHFFFAOYSA-N</cmd>

- Display a molecule by CID:
<cmd>display mol 2244</cmd>
""",
        )
    )
    statements.append(
        Forward(d_isplay + sources + (molecule_identifier)("molecule_identifier"))("display_property_sources")
    )
    grammar_help.append(
        help_dict_create(
            name="display sources",
            category="Molecules",
            command="display sources <name> | <smiles> | <inchi> | <inchikey> | <cid>",
            description=f"""
Display the sources of a molecule's properties, attributing back to how they were calculated or sourced.

{WORKING_SET_PRIORITY}

{SPECIFY_MOL}

{USING_NAME}
            

""",
        )
    )
    # ---
    # Rename molecule
    statements.append(
        Forward(
            rename
            + molecule
            + molecule_identifier("molecule_identifier")
            + a_s
            + (Word(alphas, alphanums + "_")("new_name"))
        )("rename_molecule")
    )
    grammar_help.append(
        help_dict_create(
            name="rename molecule",
            category="Molecules",
            command="rename molecule <molecule_identifer_string> as <molecule_name>",
            description="""
Rename a molecule in the current working set.

{MOL_SHORTHAND}

Example:
Let's say you've added a molecule "CC(=O)OC1=CC=CC=C1C(=O)O" to your current working set, you can then rename it as such:
<cmd>rename molecule CC(=O)OC1=CC=CC=C1C(=O)O as Aspirin</cmd>
""",
        )
    )

    # ---
    # Export molecule
    statements.append(
        Forward(
            export
            + molecule
            + (molecule_identifier)("molecule_identifier")
            + Optional(CaselessKeyword("as") + CaselessKeyword("file"))("as_file")
        )("export_molecule")
    )
    # @later - the way this works is hard to explain and not very intuitive. I would suggest making "as file" mandatory in CLI, so behavior is consistent for both CLI and Notebook.
    grammar_help.append(
        help_dict_create(
            name="export molecule",
            category="Molecules",
            command="export molecule <name> | <smiles> | <inchi> | <inchikey> | <cid> [ as file ]",
            description=f"""
When run inside a Notebook, this will return a dictionary of the molecule's properties. When run from the command line, or when `as file` is set, the molecule will be saved to your workspace as a JSON file, named after the molecule's identifying string.

{MOL_SHORTHAND}

{WORKING_SET_PRIORITY}

{USING_NAME}

Examples
- <cmd>export molecule aspirin</cmd>
- <cmd>export molecule aspirin as file</cmd>
""",
        )
    )

    # ---
    # Remove molecule
    statements.append(Forward(remove + molecule + molecule_identifier("molecule_identifier"))("remove_molecule"))
    grammar_help.append(
        help_dict_create(
            name="remove molecule",
            category="Molecules",
            command="remove molecule <name> | <smiles> | <inchi> | <inchikey> | <cid>",
            description="""
Remove a molecule from the current working set.

{MOL_SHORTHAND}
            
Examples:
- Remove a molecule by name:
<cmd>remove molecule Aspirin</cmd>

- Remove a molecule by SMILES:
<cmd>remove molecule CC(=O)OC1=CC=CC=C1C(=O)O</cmd>

- Remove a molecule by InChIKey:
<cmd>remove mol  BSYNRYMUTXBXSQ-UHFFFAOYSA-N</cmd>

- Remove a molecule by InChI:
<cmd>remove mol  InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)</cmd>

- Remove a molecule by CID:
<cmd>remove mol 2244</cmd>
""",
        )
    )

    # ---
    # List molecules
    statements.append(Forward(lister + molecules)("list_molecules"))
    grammar_help.append(
        help_dict_create(
            name="list molecules",
            category="Molecules",
            command="list molecules",
            description="List all molecules in the current working set.",
        )
    )

    statements.append(
        Forward(save + molecule_set + a_s + Word(alphas, alphanums + "_")("molecule-set_name"))("save_molecule-set")
    )

    # ---
    # Save molecules
    grammar_help.append(
        help_dict_create(
            name="save molecule-set",
            category="Molecules",
            command="save molecule-set as <molecule_set_name>",
            description="""
Save the current molecule working set to your workspace.

Example:
<cmd>save molset as my_working_set</cmd>
""",
        )
    )

    # ---
    # Load molecule set
    statements.append(
        Forward(load + molecule_set + Word(alphas, alphanums + "_")("molecule-set_name"))("load_molecule-set")
    )
    grammar_help.append(
        help_dict_create(
            name="load molecule-set",
            category="Molecules",
            command="load molecule-set|molset <molecule-set_name>",
            description="""
Load a molecule set from your workspace, and set it as your current working set.

Example:
<cmd>load molset my_working_set</cmd>
""",
        )
    )

    # ---
    # List molecule set
    statements.append(Forward(lister + molecule_sets)("list_molecule-sets"))
    grammar_help.append(
        help_dict_create(
            name="list molecule-sets",
            category="Molecules",
            command="list molecule-sets",
            description="List all molecule sets in your workspace.",
        )
    )

    # ---
    # Enrich molecule set
    statements.append(Forward(enrich + molecule_set + w_ith + analysis)("load_analysis"))
    grammar_help.append(
        help_dict_create(
            name="enrich molecule-set",
            category="Molecules",
            command="enrich molecule-set with analysis",
            description="""Enrich every molecule in your current working set with the last analysis result. This assumes that the current working set was the input or result for the analysis.

            This command currently covers results from the following Analysis commands:
            - RXN Toolkit <cmd>predict Reaction</cmd>
            - RXN Toolkit <cmd>predict retrosynthesis </cmd>
            - DS4SD Toolkit <cmd>search for patents containing molecule</cmd>
            - DS4SD Toolkit <cmd>search for similiar molecules</cmd>
              """,
        )
    )

    # ---
    # Clear analysis cache
    statements.append(Forward(clear + analysis + cache)("clear_analysis"))
    grammar_help.append(
        help_dict_create(
            name="clear analysis cache",
            category="Molecules",
            command="clear analysis cache",
            description="Clear the cache of analysis results for your current workspace.",
        )
    )
    statements.append(
        Forward(create + molecule + molecule_identifier("smiles") + name + (Word(alphas, alphanums + "_")("name")))(
            "create_molecule"
        )
    )

    # ---
    # Clear molecules
    statements.append(Forward(clear + molecules)("clear_molecules"))
    grammar_help.append(
        help_dict_create(
            name="clear Molecules",
            category="Molecules",
            command="clear molecules",
            description="Clear the working set of molecules.",
        )
    )
    statements.append(
        Forward(create + molecule + molecule_identifier("smiles") + name + (Word(alphas, alphanums + "_")("name")))(
            "create_molecule"
        )
    )

    # ---
    # Create molecule
    grammar_help.append(
        help_dict_create(
            name="create molecule",
            category="Molecules",
            command="create molecule <smiles_string> name <molecule_name>",
            description="""
Create a base molecule and add it to your working list.

Note that other identifiers (InChI and formula) will be calculated, but no other properties (like InChIKey, CID, etc.) will be fetched from PubChem.

Example:
<cmd>create molecule CC(=O)OC1=CC=CC=C1C(=O)O name my_aspirin</cmd>
""",
        )
    )

    # ---
    # Get molecule property
    statements.append(
        Forward("@" + molecule_identifier("molecule_identifier") + ">>" + mol_properties("property"))("mol_property")
    )
    grammar_help.append(
        help_dict_create(
            name="@<molecule>",
            category="Molecules",
            command="@(<name> | <smiles> | <inchi> | <inchikey> | <cid>)>><molecule_property_name>",
            description=f"""
Request a molecule's certain property.

The <cmd>@</cmd> symbol should be followed by the molecule's name, SMILES, InChI, InChIKey or CID, then after the <cmd>>></cmd> include one of the properties mentioned below.

E.g. <cmd>@aspirin>>xlogp</cmd>

{SPECIFY_MOL}

{USING_NAME}

Examples:
- Obtain the molecular weight of the molecule known as Aspirin.
<cmd>@aspirin>>molecular_weight</cmd>

- Obtain the canonical smiles string for a molecule known as Aspirin.
<cmd>@aspirin>>canonical_smiles</cmd>

- Obtain a molecules xlogp value using a SMILES string.
<cmd>@CC(=O)OC1=CC=CC=C1C(=O)O>>xlogp</cmd>

Available properties: <cmd>{'</cmd>, <cmd>'.join(m_props)}</cmd>
""",
        )
    )

    # ---
    # Load molecules from file
    statements.append(
        Forward(
            load + molecules + using + file + desc("moles_file") + Optional((merge + w_ith + pubchem))("pubchem_merge")
        )("load_molecules_file")
    )  # From mols file
    grammar_help.append(
        help_dict_create(
            name="load molecules",
            category="Molecules",
            command="load molecules using file '<csv_or_sdf_filename>' [ merge with pubchem ]",
            description="Load molecules from a CSV or SDF file into the molecule working set. Optionally you can add <cmd>merge with pubchem</cmd> to the command to fill in missing properties of the molecule.",
        )
    )

    # ---
    # Load molecules from dataframe
    statements.append(
        Forward(
            load
            + molecules
            + using
            + CaselessKeyword("dataframe")
            + molecule_identifier("in_dataframe")
            + Optional((merge + w_ith + pubchem))("pubchem_merge")
        )("load_molecules_dataframe")
    )  # From dataframe
    grammar_help.append(
        help_dict_create(
            name="load molecules",
            category="Utility",
            command="load molecules using dataframe <dataframe>",
            description="Load molecules into the molecule working set from a dataframe.",
        )
    )

    # ---
    # Export molecules
    statements.append(Forward((export + molecules + Optional(a_s + desc("csv_file_name")))("export_molecules")))
    grammar_help.append(
        help_dict_create(
            name="export molecules",
            category="Molecules",
            command="export molecules [ as <csv_filename> ]",
            description="""
Export the molecules in the current working set.

When run inside a Notebook, this will return a dataframe. When run from the command line, the molecules will be saved to your workspace as a CSV file named "result_#.csv". The rows will be numbered with the highest number representing the latest molecule that was added.
""",
        )
    )

    # ---
    # Show individual molecule in browser.
    statements.append(
        Forward(show("show") + molecule + molecule_identifier("molecule_identifier"))("show_molecule")
    )  # From mol json file
    grammar_help.append(
        help_dict_create(
            name="show molecule",
            category="Molecules",
            command="show molecule <name> | <smiles> | <inchi> | <inchikey> | <cid>",
            description="Inspect a molecule in the browser.",
        )
    )

    # ---
    # Show molecule-set in browser.
    statements.append(
        Forward(
            show("show")
            + molecules
            + using
            + file
            # + desc("moles_file")  # From mols file
            + Word(alphas, alphanums + "_")("molecule-set_name")
            + Optional(a_s + CaselessKeyword("molsobject")("object"))  # Return as molsobject
            + Optional(save + a_s + desc("results_file"))  # Save as csv/sdf
        )("show_molecules")
    )
    grammar_help.append(
        help_dict_create(
            name="show molecules",
            category="Molecules",
            command="show molecules using ( file '<mols_file>' | dataframe <dataframe> ) [ save as '<sdf_or_csv_file>' | as molsobject ]",
            description=f"""Launch the molecule viewer { 'in your browser ' if is_notebook_mode() else '' }to examine and select molecules from a SMILES sdf/csv dataset.

    Examples:
    - <cmd>show molecules using file 'base_molecules.sdf' as molsobject</cmd>
    - <cmd>show molecules using dataframe my_dataframe save as 'selection.sdf'</cmd>
    """,
        )
    )
