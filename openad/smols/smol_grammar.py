from openad.core.help import help_dict_create
from openad.smols.smol_functions import SMOL_PROPERTIES as m_props
from openad.helpers.general import is_notebook_mode

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
_mol_properties.extend(m_props)
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
SUPPORTED_FILE_FORMATS = """Supported file formats:
- molset (.molset.json)
- SDF (.sdf)
- CSV (.csv)
- SMILES (.smi)"""
SPECIFY_MOL = "You can specify any molecule by SMILES or InChI, and PubChem classified molecules also by name, InChIKey or their PubChem CID. \n A molecule identifier can be in single quotes or defined with unquoted text. If you have spaces in your molecule identifier e.g. a name, then you must user a single quoted string"
WORKING_SET_PRIORITY = "If the requested molecule exists in your current working set, that version will be used."
USING_NAME = "If you use the name of a molecule, the tool will do a caseless search of the names and synonyms first in current molecule working set, then on PubChem."


def smol_grammar_add(statements, grammar_help):
    """defines the grammar available for managing molecules"""

    # ---
    # Add molecule
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
            command="add mol|molecule <name> | <smiles> | <inchi> | <inchikey> | <cid> [ as '<name>' ] [ basic ] [ force ]",
            description=f"""
This command is how you add a molecule to your current molecule working set in memory. When adding a molecule by name, this name will become the molecule's identifying string. 

It will take any molecules identifier from the following categories:
    -<cmd>smiles </cmd>
    -<cmd>name or synonym</cmd>
    -<cmd>smiles</cmd>
    -<cmd>inchi</cmd>
    -<cmd>inchikey </cmd>
    -<cmd>cid </cmd>

Options :
    - <cmd>as <name> </cmd>: if the <cmd> as '<name>' </cmd> not used the molecule the  molecule identfier will be used for the molecules name. if the <cmd> as '<name>' </cmd> not used the molecule the  molecule identfier will be used for the molecules name.
        You can set or override an name later for  any molecule with the <cmd>rename molecule</cmd> command.
    - <cmd> basic </cmd> Creates a molecule that does not have its properties and synonyms populated with pubchem data, this feature is only valid with a SMILES molecule identifier
    - <cmd>force</cmd>: The <cmd>force</cmd> option allows you to ovveride the confirmation that you wish to add a molecule.



    
{MOL_SHORTHAND}

{SPECIFY_MOL}

{USING_NAME}


Examples of how to add a molecule to your molecule working set:
- Add a molecule by name:
<cmd>add molecule aspirin</cmd>
or with single quotes
<cmd> display molecule 'Aspirin 325 mg' </cmd>

- Add a molecule by name and force through the acknowledgement to add it:
<cmd>add molecule aspirin force</cmd>

- Add a molecule by SMILES:
<cmd>add molecule CC(=O)OC1=CC=CC=C1C(=O)O</cmd>

- Add a molecule by SMILES without populated pubchem properties:
<cmd>add molecule CC(=O)OC1=CC=CC=C1C(=O)O basic </cmd>

- Add a molecule by CID:
<cmd>add mol 2244</cmd>

- Add a molecule by InChI:
<cmd>add mol InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)</cmd>
 or with single quotes
 <cmd>add mol 'InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)'</cmd>

- Add a molecule by InChIKey:
<cmd>add mol BSYNRYMUTXBXSQ-UHFFFAOYSA-N</cmd>

- Add a molecule by InChIKey nd set its name to "mymol":
<cmd>add mol BSYNRYMUTXBXSQ-UHFFFAOYSA-N as 'mymol' </cmd>

- Add a molecule by SMILES nd set its name to "mymol" and not prepopulate values from pubchem:
<cmd>add mol CC(=O)OC1=CC=CC=C1C(=O)O as 'mymol' basic </cmd>
""",
        )
    )

    # ---
    # Display molecule
    statements.append(
        Forward(d_isplay + molecule + (molecule_identifier | desc)("molecule_identifier"))("display_molecule")
    )
    grammar_help.append(
        help_dict_create(
            name="display molecule",
            category="Small Molecules",
            command="display molecule <name> | <smiles> | <inchi> | <inchikey> | <cid>",
            description=f"""
This command will display a molecule's identifiers, propoerties, synonyms and any Analysis results it has been enriched with.
if the molecule exists in the current molecule workling list in memory the molecule will be retrieved from there if not pubchem will be checked to see if the molecule and its information is avialable there.

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
        Forward(d_isplay + sources + (molecule_identifier | desc)("molecule_identifier"))("display_property_sources")
    )
    grammar_help.append(
        help_dict_create(
            name="display sources",
            category="Molecule Working Set",
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
            + (molecule_identifier | desc)("molecule_identifier")
            + a_s
            + (Word(alphas, alphanums + "_")("new_name"))
        )("rename_molecule")
    )
    grammar_help.append(
        help_dict_create(
            name="rename molecule",
            category="Molecule Working Set",
            command="rename molecule <molecule_identifer_string> as <molecule_name>",
            description=f"""
Rename a molecule in the current working set.

{MOL_SHORTHAND}

Example:
Let's say you've added a molecule "CC(=O)OC1=CC=CC=C1C(=O)O" to your current molecule working set, you can then rename it as such:
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
            + (molecule_identifier | desc)("molecule_identifier")
            + Optional(CaselessKeyword("as") + CaselessKeyword("file"))("as_file")
        )("export_molecule")
    )
    # @later - the way this works is hard to explain and not very intuitive. I would suggest making "as file" mandatory in CLI, so behavior is consistent for both CLI and Notebook.
    grammar_help.append(
        help_dict_create(
            name="export molecule",
            category="Molecule Working Set",
            command="export molecule <name> | <smiles> | <inchi> | <inchikey> | <cid> [ as file ]",
            description=f"""
When run inside a jupyter lab notebook, this will return a dictionary of the molecule's properties. When run from the command line, or when <cmd>as file</cmd> is set, the molecule will be saved to your workspace as a JSON file, named after the molecule's identifying string.
If the molecule is in your current working set it will be retrieved from there, if the molecule is not there pubchem will be called to retrieve the molecule.

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
    statements.append(
        Forward(remove + molecule + (molecule_identifier | desc)("molecule_identifier") + Optional(force)("force"))(
            "remove_molecule"
        )
    )
    grammar_help.append(
        help_dict_create(
            name="remove molecule",
            category="Molecule Working Set",
            command="remove mol|molecule <name> | <smiles> | <inchi> | <inchikey> | <cid> [ force ]",
            description=f"""
Remove a molecule from the current working set based on a given molecule identifier.

{MOL_SHORTHAND}
            
Examples:
- Remove a molecule by name:
<cmd>remove molecule Aspirin</cmd>

- Remove a molecule by SMILES:
<cmd>remove molecule CC(=O)OC1=CC=CC=C1C(=O)O</cmd>

- Remove a molecule by InChIKey:
<cmd>remove mol  BSYNRYMUTXBXSQ-UHFFFAOYSA-N</cmd>

- Remove a molecule by InChI:
<cmd>remove mol  InChI='1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)'</cmd>

- Remove a molecule by CID:
<cmd>remove mol 2244</cmd>
""",
        )
    )

    # ---
    # List molecules
    statements.append(Forward(l_ist + molecules)("list_molecules"))
    grammar_help.append(
        help_dict_create(
            name="list molecules",
            category="Molecule Working Set",
            command="list molecules",
            description="List all molecules in the current working set.",
        )
    )

    # ---
    # Display my-mols (working set) in browser/iFrame (CLI/Jupyter)
    statements.append(Forward(show("show") + molecules)("show_molecules"))
    grammar_help.append(
        help_dict_create(
            name="show molecules",
            category="Molecule Working Set",
            command="show mols|molecules",
            description="Display the current molecule working set in the GUI.",
        )
    )

    # ---
    # Save molecules
    statements.append(
        Forward(save + molecule_set + a_s + Word(alphas, alphanums + "_")("molecule-set_name"))("save_molecule-set")
    )
    grammar_help.append(
        help_dict_create(
            name="save molecule-set",
            category="Molecule Sets",
            command="save molecule-set as <molecule_set_name>",
            description="""
Save the current molecule workking list to a molecule-set in your workspace.

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
            category="Molecule Sets",
            command="load molecule-set|molset <molecule-set_name>",
            description="""
Loads a molecule-set from your workspace, and replaces your current list of molecules with the molecules from the given  molecule-set.
Example:
<cmd>load molset my_working_set</cmd>
""",
        )
    )

    # Load molecule set
    statements.append(
        Forward(
            merge
            + molecule_set
            + Word(alphas, alphanums + "_")("molecule-set_name")
            + Optional(merge + only)("merge_only")
            + Optional(append + only)("append_only")
        )("merge_molecule-set")
    )
    grammar_help.append(
        help_dict_create(
            name="merge molecule-set",
            category="Molecule Sets",
            command="merge molecule-set|molset <molecule-set_name> [merge only] [append only]",
            description="""
This command merges a molecule-set from your workspace into cour current molecule working set in memory, and updates properties/Analysis in existing molecules or appends new molecules to the working set.

Options:
    - <cmd> merge only</cmd> Only merges with existing molecules in list
    - <cmd> append only</cmd> Only append molecules not in list
<cmd>merge molset my_working_set</cmd>

""",
        )
    )

    # ---
    # List molecule set
    statements.append(Forward(l_ist + molecule_sets)("list_molecule-sets"))
    grammar_help.append(
        help_dict_create(
            name="list molecule-sets",
            category="Molecule Sets",
            command="list molecule-sets",
            description="List all molecule sets in your workspace.",
        )
    )

    # ---
    # Enrich molecule set
    statements.append(Forward(enrich + molecules + w_ith + analysis)("load_analysis"))
    grammar_help.append(
        help_dict_create(
            name="enrich molecules",
            category="Molecule Working Set",
            command="enrich molecules with analysis",
            description="""
Enrich the molecules in your working set with the results of the last performed analysis.
This assumes that your molecule working set contains either the input molecule or any of the result molecules from the analysis.

Currently supported analysis commands:

RXN:
- <cmd>predict reaction</cmd>
- <cmd>predict retrosynthesis </cmd>

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
            description="this command clears the cache of analysis results for your current workspace.",
        )
    )

    # ---
    # Clear molecules
    statements.append(Forward(clear + molecules)("clear_molecules"))
    grammar_help.append(
        help_dict_create(
            name="clear Molecules",
            category="Molecule Working Set",
            command="clear molecules",
            description="This command clears the molecule working set.",
        )
    )

    # ---
    # Get molecule property
    statements.append(
        Forward("@" + (molecule_identifier | desc)("molecule_identifier") + ">>" + mol_properties("property"))(
            "get_smol_prop"
        )
    )
    grammar_help.append(
        help_dict_create(
            name="@<molecule>",
            category="Small Molecules",
            command="@(<name> | <smiles> | <inchi> | <inchikey> | <cid>)>><molecule_property_name>",
            description=f"""This command request the given property of a molecule, it will first try and retrieve the provided molecule from your molecule working set, if it is not there it will will try and retrieve the molecule from pubchem.

The <cmd>@</cmd> symbol should be followed by the molecule's name, SMILES, InChI, InChIKey or CID, then after the <cmd>>></cmd> include one of the properties mentioned below.

E.g. <cmd>@aspirin>>xlogp</cmd>

{SPECIFY_MOL}

{USING_NAME}

Examples of how to retrieve the value of a molecules property:
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
            load
            + molecules
            + f_rom
            + file
            + desc("moles_file")
            + Optional((merge + w_ith + pubchem))("pubchem_merge")
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
            + Optional((merge + w_ith + pubchem))("pubchem_merge")
            + Optional(CaselessKeyword("append"))("append")
        )("load_molecules_file-DEPRECATED")
    )
    grammar_help.append(
        help_dict_create(
            name="load molecules",
            category="Molecule Working Set",
            command="load molecules from file '<filename.molset.json|sdf|csv|smi>' [ merge with pubchem ] [ append ]",
            description=f"""Load molecules from a file into your molecule working set.

{MOLS_SHORTHAND}

{SUPPORTED_FILE_FORMATS}

Options:
- Append <cmd>merge with pubchem</cmd> to enrich the molecule with data from pubchem
- Append <cmd>append</cmd> to append the molecules to the existing working set instead of overwriting it

Examples:
- <cmd>load molecules from file 'my_molecules.molset.json'</cmd>
- <cmd>load mols from file 'my_molecules.sdf'` appen</cmd>
- <cmd>load molecules from file 'my_molecules.csv'</cmd>
- <cmd>load mols from file 'my_molecules.smi'` merge with pubche</cmd>
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
            + Optional((merge + w_ith + pubchem))("pubchem_merge")
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
            + Optional((merge + w_ith + pubchem))("pubchem_merge")
            + Optional(CaselessKeyword("append"))("append")
        )("load_molecules_dataframe-DEPRECATED")
    )
    grammar_help.append(
        help_dict_create(
            name="load Molecule Working Set",
            category="Molecule Working Set",
            command="load molecules from dataframe <dataframe> [ merge with pubchem ] [ append ]",
            description=f"""Load molecules from a dataframe into your molecule working set.

{MOLS_SHORTHAND}

This command only works when called from a Jupyter Notebook or the API.

Options:
- Append <cmd>merge with pubchem</cmd> to enrich the molecule with data from pubchem
- Append <cmd>append</cmd> to append the molecules to the existing working set instead of overwriting it

Examples:
- <cmd>load molecules from dataframe my_dataframe</cmd>
- <cmd>load mols from dataframe my_dataframe append</cmd>
- <cmd>load mols from dataframe my_dataframe merge with pubchem</cmd>
""",
        )
    )

    statements.append(
        Forward(
            CaselessKeyword("merge")
            + molecules
            + data
            + using
            + CaselessKeyword("dataframe")
            + molecule_identifier("in_dataframe")
            + Optional((merge + w_ith + pubchem))("pubchem_merge")
        )("merge_molecules_data_dataframe")
    )  # From dataframe
    grammar_help.append(
        help_dict_create(
            name="merge molecules data",
            category="Utility",
            command="merge molecules data using dataframe <dataframe> [ merge with pubchem ]",
            description="""            

            
This command merges molecules into the molecule working set from a dataframe. 
     
It takes files with the columns named: 

<cmd>subject or <cmd>smiles</cmd>: molecules similes string
<cmd>property</cmd> : the property generation name
<cmd>result</cmd> : the value of the property

Sample input file

subject                                                               property                        result
--------------------------------------------------------------------  -------------------------  -----------
O=C(N)C(F)(OC(F)(F)C(F)(F)C(F)(F)F)C(F)(F)F                           molecular_weight               329.065
O=C(O)C(F)(NC(F)(F)C(F)(F)C(F)(F)F)C(F)(F)F                           molecular_weight               329.065
O=C(O)C(F)(OC(F)(F)C(F)(F)C)CF                                        molecular_weight               240.099
O=C(O)C(F)(OC(O)(F)C(F)(F)C(F)(F)F)C(F)(F)F                           molecular_weight               328.058
O=C(O)C(F)(OC(Cl)(F)C(F)(F)C(F)(F)F)C(F)(F)F                          molecular_weight               346.504
O=C(O)C(F)(OC(F)(F)C(F)(O)C(F)(F)F)C(F)(F)F                           molecular_weight               328.058
O=C(O)C(F)(OC(F)(O)C(F)(F)C(F)(F)F)C(F)(F)F                           molecular_weight               328.058
O=C(O)C(F)(OC(F)(F)C(F)(Br)C(F)(F)F)C(F)(F)F                          molecular_weight               390.955
O=C(O)C(F)OC(O)(F)C(F)(F)C(F)(F)F                                     molecular_weight               260.061


Example Command:

merge molecules from a data frame called <cmd>new_props</cmd>

- <cmd> merge molecules data using dataframe new_props</cmd>

to perform the same load and merge with pubchem data simply add the <cmd> merge with pubchem </cmd> clause to the end of the command 

- <cmd> merge molecules data using dataframe new_props merge with pubchem</cmd> """,
        )
    )

    # ---
    # Export molecules
    statements.append(Forward((export + molecules + Optional(a_s + desc("csv_file_name")))("export_mws")))
    grammar_help.append(
        help_dict_create(
            name="export molecules",
            category="Molecule Working Set",
            command="export molecules [ as '<filename.molset.json|sdf|csv|smi>' ]",
            description=f"""Export your molecule working set as a dataframe (Jupyter/API) or a file (CLI).

{MOLS_SHORTHAND}

When run inside a Jupyter Notebook or from the API, the <cmd>as <filename></cmd> clause will be ignored and a dataframe will be returned.

When exporting as a file, the filename's extension will define what format the molecule are exported as.

{SUPPORTED_FILE_FORMATS}

If no filename or extension is provided, the molecules will be saved as molset file.
""",
        )
    )

    # ---
    # Show individual molecule in browser.
    statements.append(
        Forward(show("show") + molecule + (molecule_identifier | desc)("molecule_identifier"))("show_mol")
    )  # From mol json file
    grammar_help.append(
        help_dict_create(
            name="show molecule",
            category="Small Molecules",
            command="show mol|molecule <name> | <smiles> | <inchi> | <inchikey> | <cid>",
            description=f"""
Inspect a molecule in the browser. If a molecule is not in the current Molecule Working set it will pull the result from Pubchem.

{MOL_SHORTHAND}

When you show a molecule by SMILES or InChI, we can display it immediately. When you show a molecule by name, InChIKey or PubChem CID, we need to first retrieve it from PubChem, which can take a few seconds.

Examples:
- <cmd>show mol aspirin</cmd>
- <cmd>show mol CC(=O)OC1=CC=CC=C1C(=O)O</cmd>
- <cmd>show mol InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)</cmd>
- <cmd>show mol 2244</cmd>
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
            category="Molecule Sets",
            command="show molset|molecule-set '<molset_or_sdf_or_smi_path>' | using dataframe <dataframe>",
            description=f"""
Launch the molset viewer { 'in your browser ' if is_notebook_mode() else '' }to visualize your molecule dataset.

Examples:
- <cmd>show molset 'neurotransmitters.mol.json'</cmd>
- <cmd>show molset 'neurotransmitters.sdf'</cmd>
- <cmd>show molset 'neurotransmitters.smi'</cmd>
- <cmd>show molset my_dataframe</cmd>
""",
        )
    )
