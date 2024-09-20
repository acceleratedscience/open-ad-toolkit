from openad.core.help import help_dict_create
from openad.smols.smol_functions import MOL_PROPERTIES as m_props
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
MOLSET_SHORTHAND = "You can use the 'molset' shorthand instead of 'molecule-set'."
SPECIFY_MOL = "You can specify any molecule by SMILES or InChI, and PubChem classified molecules also by name, InChIKey or their PubChem CID. \n A molecule identifier can be in single quotes or defined with unquoted text. If you have spaces in your molecule identifier e.g. a name, then you must user a single quoted string"
WORKING_SET_PRIORITY = "If the requested molecule exists in your current working list, that version will be used."
USING_NAME = "If you use the name of a molecule, the tool will do a caseless search of the names and synonyms first in current working list, then on pubchem."


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
            category="Molecules Working Set",
            command="add mol|molecule <name> | <smiles> | <inchi> | <inchikey> | <cid> [ as '<name>' ] [ basic ] [ force ]",
            description=f"""
This command is how you add a molecule to a current working list of molecules in memory. When adding a molecule by name, this name will become the molecule's identifying string. 

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


Examples of how to add a molecule to your working list:
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
            category="Molecules Working Set",
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
            category="Molecules Working Set",
            command="rename molecule <molecule_identifer_string> as <molecule_name>",
            description=f"""
Rename a molecule in the current working list.

{MOL_SHORTHAND}

Example:
Let's say you've added a molecule "CC(=O)OC1=CC=CC=C1C(=O)O" to your current working list of molecules, you can then rename it as such:
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
            category="Molecules Working Set",
            command="export molecule <name> | <smiles> | <inchi> | <inchikey> | <cid> [ as file ]",
            description=f"""
When run inside a jupyter lab notebook, this will return a dictionary of the molecule's properties. When run from the command line, or when <cmd>as file</cmd> is set, the molecule will be saved to your workspace as a JSON file, named after the molecule's identifying string.
If the molecule is in your current working list it will be retrieved from there, if the molecule is not there pubchem will be called to retrieve the molecule.

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
            category="Molecules Working Set",
            command="remove mol|molecule <name> | <smiles> | <inchi> | <inchikey> | <cid> [ force ]",
            description=f"""
Remove a molecule from the current working list based on a given molecule identifier.

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
            category="Molecules Working Set",
            command="list molecules",
            description="List all molecules in the current working list.",
        )
    )

    # ---
    # Display my-mols (working set) in browser/iFrame (CLI/Jupyter)
    statements.append(Forward(show("show") + molecules)("show_molecules"))
    grammar_help.append(
        help_dict_create(
            name="show molecules",
            category="Molecules Working Set",
            command="show mols|molecules",
            description="Display the current working list of molecules in the GUI.",
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
This command merges a molecule-set from your workspace into cour current working list of molecules in memory, and updates properties/Analysis in existing molecules or appends new molecules to the working list.

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
            category="Molecules Working Set",
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
            category="Molecules Working Set",
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
            category="Molecules Working Set",
            command="clear molecules",
            description="This command clears the working list of molecules.",
        )
    )

    # ---
    # Get molecule property
    statements.append(
        Forward("@" + (molecule_identifier | desc)("molecule_identifier") + ">>" + mol_properties("property"))(
            "mol_property"
        )
    )
    grammar_help.append(
        help_dict_create(
            name="@<molecule>",
            category="Small Molecules",
            command="@(<name> | <smiles> | <inchi> | <inchikey> | <cid>)>><molecule_property_name>",
            description=f"""
This command request the given property of a molecule, it will first try and retrieve the provided molecule from your working list of molecules, if it is not there it will will try and retrieve the molecule from pubchem.

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
            + using
            + file
            + desc("moles_file")
            + Optional((merge + w_ith + pubchem))("pubchem_merge")
            + Optional(CaselessKeyword("append"))("append")
        )("load_molecules_file")
    )  # From mols file
    grammar_help.append(
        help_dict_create(
            name="load molecules",
            category="Molecules Working Set",
            command="load molecules using file '<csv_or_sdf_filename>' [ merge with pubchem ] [append]",
            description="""This command Loads molecules from a CSV or SDF file into the molecule working list. 
            
            Options:
             - you can add <cmd>merge with pubchem</cmd> to the command to fill in missing properties of the molecule.
             - you can append to the existing working set using the command <cmd> append </append> 

             and example data set is as follows

    cid      SMILES                                                                 chemical_name                    molecular weight    xlogp3    
--------   ---------------------------------------------------------------------  --------------------------------  ------------------  --------  
  114481  C(=O)(C(C(F)(F)F)(OC(C(C(F)(F)F)(F)F)(F)F)F)O                          propanoic acid                     330.05              3.6    
   67821  C(=O)(C(C(C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)O      Perfluorononanoic acid             464.08              5.6     
    9554  C(=O)(C(C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)O             Perfluorooctanoic acid             414.07              4.9    
   74483  C(C(C(C(C(F)(F)S(=O)(=O)O)(F)F)(F)F)(F)F)(C(C(C(F)(F)F)(F)F)(F)F)(F)F  Perfluorooctane sulfonic acid      500.13               5     
   67734  C(C(C(C(F)(F)S(=O)(=O)O)(F)F)(F)F)(C(C(F)(F)F)(F)F)(F)F                Perfluorohexanesulphonic acid      400.12              3.7     
16760155  C(C(C(C(F)(F)S(=O)(=O)[O-])(F)F)(F)F)(C(C(F)(F)F)(F)F)(F)F                                                399.11              3.6 
             
             
             """,
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
            + Optional(CaselessKeyword("append"))("append")
        )("load_molecules_dataframe")
    )
    # From dataframe
    grammar_help.append(
        help_dict_create(
            name="load Molecules Working Set",
            category="Molecules Working Set",
            command="load molecules using dataframe <dataframe> [ merge with pubchem ] [append]",
            description=""""            
This command Load molecules into the molecule working list from a dataframe. 

            Options:
             - you can add <cmd>merge with pubchem</cmd> to the command to fill in missing properties of the molecule. NOTE:  this will slow the process down
             - you can append to the existing working set using the command <cmd> append </append> 
            
            an example data set that is compatible is as follows


    cid      SMILES                                                                 chemical_name                    molecular weight    xlogp3    
--------   ---------------------------------------------------------------------  --------------------------------  ------------------  --------  
  114481  C(=O)(C(C(F)(F)F)(OC(C(C(F)(F)F)(F)F)(F)F)F)O                          propanoic acid                     330.05              3.6    
   67821  C(=O)(C(C(C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)O      Perfluorononanoic acid             464.08              5.6     
    9554  C(=O)(C(C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)O             Perfluorooctanoic acid             414.07              4.9    
   74483  C(C(C(C(C(F)(F)S(=O)(=O)O)(F)F)(F)F)(F)F)(C(C(C(F)(F)F)(F)F)(F)F)(F)F  Perfluorooctane sulfonic acid      500.13               5     
   67734  C(C(C(C(F)(F)S(=O)(=O)O)(F)F)(F)F)(C(C(F)(F)F)(F)F)(F)F                Perfluorohexanesulphonic acid      400.12              3.7     
16760155  C(C(C(C(F)(F)S(=O)(=O)[O-])(F)F)(F)F)(C(C(F)(F)F)(F)F)(F)F                                                399.11              3.6 
             
             
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

            
This command merges molecules into the molecule working list from a dataframe. 
     
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
    statements.append(Forward((export + molecules + Optional(a_s + desc("csv_file_name")))("export_molecules")))
    grammar_help.append(
        help_dict_create(
            name="export molecules",
            category="Molecules Working Set",
            command="export molecules [ as <csv_filename> ]",
            description="""
This command exports the molecules in the current working list of molecules.

When run inside a Notebook, this will return a dataframe. When run from the command line, the molecules will be saved to your workspace as a CSV file named "result_#.csv". The rows will be numbered with the highest number representing the latest molecule that was added.
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

    # --- TRASH / DEPRECATED
    # Show molset in browser
    # ONLY HERE TO COMPARE OLD MOLS2GRID TO NEW GUI MOLSET VIEWER
    # Note: we don't allow dashes in dataframe names because it's a substraction operator and causes issues in Jupyter.
    statements.append(
        Forward(
            show("show")
            + molecules
            + using
            + CaselessKeyword("dataframe")
            + Word(alphas, alphanums + "_")("in_dataframe")  # From dataframe
            + Optional(a_s + CaselessKeyword("molsobject")("object"))  # Return as molsobject
            + Optional(save + a_s + desc("results_file"))  # Save as csv/sdf
        )("show_molsgrid_df")
    )

    # --- TRASH / DEPRECATED
    # Show molset in browser
    # ONLY HERE TO COMPARE OLD MOLS2GRID TO NEW GUI MOLSET VIEWER
    statements.append(
        Forward(
            show("show")
            + molecules
            + using
            + file
            + desc("moles_file")  # From mols file
            # + Word(alphas, alphanums + "_")("molecule-set_name")
            + Optional(a_s + CaselessKeyword("molsobject")("object"))  # Return as molsobject
            + Optional(save + a_s + desc("results_file"))  # Save as csv/sdf
        )("show_molsgrid")
    )
    # Removed from the help display
    # grammar_help.append(
    #     help_dict_create(
    #         name="show molecules",
    #         category="Small Molecules",
    #         command="show molecules using ( file '<mols_file>' | dataframe <dataframe> ) [ save as '<sdf_or_csv_file>' | as molsobject ]",
    #         description=f"""Launch the molecule viewer { 'in your browser ' if is_notebook_mode() else '' }to examine and select molecules from a SMILES sdf/csv dataset.
    # { 'if you are working from a notebook, the <cmd> as mols object </cmd> clause allows you to display the mols2grid object and use the <cmd>.get_selection()</cmd>  method to retrieve selected molecules ' if is_notebook_mode() else '' }
    # Examples of how to show molecules in mols2grid:
    # - <cmd>show molecules using file 'base_molecules.sdf' as molsobject</cmd>
    # - <cmd>show molecules using dataframe my_dataframe save as 'selection.sdf'</cmd>
    # """,
    #     )
    # )
