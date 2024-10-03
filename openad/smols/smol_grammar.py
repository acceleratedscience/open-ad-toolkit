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
DELETE_____SPECIFY_MOL = "You can specify any molecule by SMILES or InChI, and PubChem classified molecules also by name, InChIKey or their PubChem CID. \n A molecule identifier can be in single quotes or defined with unquoted text. If you have spaces in your molecule identifier e.g. a name, then you must user a single quoted string"
DELETE_____USING_NAME = "If you use the name of a molecule, the tool will do a caseless search of the names and synonyms first in current molecule working set, then on PubChem."


def smol_grammar_add(statements, grammar_help):
    """
    Grammar for managing small molecules.
    """

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
            command="add mol|molecule <name> | <smiles> | <inchi> | <inchikey> | <cid> [ as <name> ] [ basic ] [ force ]",
            description=f"""Add a molecule to your current molecule working set.

{SUPPORTED_IDENTIFIERS_BASIC}

Options:
- <cmd>as <name></cmd>: Provide a custom name for the molecule, which will be used by the software whenever refering to it going forward.
  Note: you can always update a molecule's name later by running <cmd>rename molecule <name></cmd>.
- <cmd>basic</cmd>: Create a minimal molecule without enriching it with PubChem data. This is only relevant when using a SMILES or InChI string as identifier. Because no API calls are made, this is much faster than the default behavior.
- <cmd>force</cmd>: This suppressed the confirmation step after adding a molecule, which may be desired in batch operations.

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
            command="remove mol|molecule <name> | <smiles> | <inchi> | <inchikey> | <cid> [ force ]",
            description=f"""Remove a molecule from the current working set based on a given molecule identifier.

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
            command="list mols|molecules",
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
            command="show mols|molecules",
            description=f"""Visualize the current molecule working set in an iframe (Jupyter Notebook) or in the browser (CLI).

Notes:
- {MOLS_SHORTHAND}
""",
        )
    )

    # ---
    # Enrich molecule set
    statements.append(Forward(enrich + molecules + w_ith + analysis)("load_analysis"))
    grammar_help.append(
        help_dict_create(
            name="enrich molecules",
            category="Molecule Working Set",
            command="enrich mols|molecules with analysis",
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

Please refer to the <cmd>enrich mols|molecules with analysis</cmd> command for more information about analysis results.
""",
        )
    )

    # ---
    # Display working set molecule's sources
    statements.append(
        Forward(d_isplay + molecule + (molecule_identifier | desc)("molecule_identifier"))("display_molecule")
    )
    grammar_help.append(
        help_dict_create(
            name="display molecule",
            category="Small Molecules",
            command="display mol|molecule <name> | <smiles> | <inchi> | <inchikey> | <cid>",
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
    statements.append(
        Forward(d_isplay + sources + (molecule_identifier | desc)("molecule_identifier"))("display_property_sources")
    )
    grammar_help.append(
        help_dict_create(
            name="display sources",
            category="Molecule Working Set",
            command="display sources <name> | <smiles> | <inchi> | <inchikey> | <cid>",
            description=f"""Display the sources of a molecule's properties, attributing how they were calculated or where they were sourced.

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
            command="rename mol|molecule <molecule_identifer_string> as <molecule_name>",
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
            + Optional(enrich)("pubchem_merge")
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
            + Optional((merge + w_ith + pubchem))("pubchem_merge")  # <-- changed
            + Optional(CaselessKeyword("append"))("append")
        )("load_molecules_file-DEPRECATED")
    )
    grammar_help.append(
        help_dict_create(
            name="load molecules",
            category="Molecule Working Set",
            command="load mols|molecules from file '<filename.molset.json|sdf|csv|smi>' [ enrich ] [ append ]",
            description=f"""Load molecules from a file into your molecule working set.

{SUPPORTED_FILE_FORMATS}

Options:
- Append <cmd>enrich</cmd> to enrich the molecule with data from pubchem.
- Append <cmd>append</cmd> to append the molecules to the existing working set instead of overwriting it.

Notes:
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
            + Optional(enrich)("pubchem_merge")
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
            command="load mols|molecules from dataframe <dataframe> [ enrich ] [ append ]",
            description=f"""Load molecules from a dataframe into your molecule working set.

Options:
- Append <cmd>enrich</cmd> to enrich the molecule with data from pubchem
- Append <cmd>append</cmd> to append the molecules to the existing working set instead of overwriting it

Notes:
- This command only works when called from a Jupyter Notebook or the API.
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
            + Optional(enrich)("pubchem_merge")
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
            + Optional((merge + w_ith + pubchem))("pubchem_merge")  # <-- changed
        )("merge_molecules_data_dataframe-DEPRECATED")
    )
    grammar_help.append(
        help_dict_create(
            name="merge molecules data",
            category="Molecule Working Set",
            command="merge mols|molecules data from dataframe <dataframe> [ enrich ]",
            description="""Merges molecule data from a dataframe into the molecules in your working set.
    
It takes files with columns named as follows:
- <cmd>subject</cmd> or <cmd>smiles</cmd>: molecules similes string
- <cmd>property</cmd>: the name of the property to be merged
- <cmd>result</cmd>: the value of the property to be nmerged

Sample input file:

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


Examples:

- Merge molecule data from a dataframe called <cmd>new_props</cmd>:
  <cmd>merge molecules data from dataframe new_props</cmd>

- Merge molecule data from a dataframe called <cmd>new_props</cmd>, while enriching the molecules with PubChem data:
  <cmd>merge molecules data from dataframe new_props enrich</cmd>
""",
        )
    )

    # ---
    # Export single molecule from working set
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
            command="export mol|molecule <name> | <smiles> | <inchi> | <inchikey> | <cid> [ as file ]",
            description=f"""
When run inside a jupyter lab notebook, this will return a dictionary of the molecule's properties. When run from the command line, or when <cmd>as file</cmd> is set, the molecule will be saved to your workspace as a JSON file, named after the molecule's identifying string.
If the molecule is in your current working set it will be retrieved from there, if the molecule is not there pubchem will be called to retrieve the molecule.

{MOL_SHORTHAND}

{MOL_LOOKUP_PRIORITY}

{DELETE_____USING_NAME}

Examples
- <cmd>export molecule aspirin</cmd>
- <cmd>export mol aspirin as file</cmd>
""",
        )
    )

    # ---
    # Export all molecules ferom working set
    statements.append(Forward((export + molecules + Optional(a_s + desc("csv_file_name")))("export_mws")))
    grammar_help.append(
        help_dict_create(
            name="export molecules",
            category="Molecule Working Set",
            command="export mols|molecules [ as '<filename.molset.json|sdf|csv|smi>' ]",
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
    # Clear molecule working set
    statements.append(Forward(clear + molecules)("clear_molecules"))
    grammar_help.append(
        help_dict_create(
            name="clear Molecules",
            category="Molecule Working Set",
            command="clear mols|molecules",
            description="This command clears the molecule working set.",
        )
    )

    #
    #
    # Molecule Sets
    #
    #

    # ---
    # Save molecules
    statements.append(
        Forward(save + molecule_set + a_s + Word(alphas, alphanums + "_")("molset_name"))("save_molecule-set")
    )
    grammar_help.append(
        help_dict_create(
            name="save molecule-set",
            category="Molecule Sets",
            command="save molset|molecule-set as <molset_name>",
            description=f"""
Save the current molecule working set to a molecule-set in your workspace.

Notes:
- {MOLSET_SHORTHAND}

Example:
<cmd>save molset as my_working_set</cmd>
""",
        )
    )

    # ---
    # Load molecule set
    statements.append(Forward(load + molecule_set + Word(alphas, alphanums + "_")("molset_name"))("load_molecule-set"))
    grammar_help.append(
        help_dict_create(
            name="load molecule-set",
            category="Molecule Sets",
            command="load molset|molecule-set <molset_name>",
            description="""
Load a molecule-set from your workspace into your working set, replacing your current list of molecules.
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
            + Word(alphas, alphanums + "_")("molset_name")
            + Optional(merge + only)("merge_only")
            + Optional(append + only)("append_only")
        )("merge_molecule-set")
    )
    grammar_help.append(
        help_dict_create(
            name="merge molecule-set",
            category="Molecule Sets",
            command="merge molset|molecule-set <molset_name> [merge only] [append only]",
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
            command="list molsets|molecule-sets",
            description="List all molecule sets in your workspace.",
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

    #
    #
    # Small Molecules
    #
    #

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

{DELETE_____SPECIFY_MOL}

{DELETE_____USING_NAME}

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
