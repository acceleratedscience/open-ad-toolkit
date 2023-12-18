from openad.core.help import help_dict_create
from openad.molecules.mol_functions import MOL_PROPERTIES as m_props
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
) = map(
    CaselessKeyword,
    "get list description using create set unset workspace workspaces context jobs exec\
          as optimize with toolkits toolkit gpu experiment add run save runs show \
              file display history data remove result from inchi inchikey smiles formula name last load results export create rename merge pubchem".split(),
)
mol = ["molecule", "mol"]
mols = ["molecules", "mols"]
molset = ["molecule-set", "molset"]
molsets = ["molecule-sets", "molsets"]
clear = CaselessKeyword("clear")
cache = CaselessKeyword("cache")
analysis = CaselessKeyword("analysis")
enrich = CaselessKeyword("enrich")
mol_properties = MatchFirst(map(CaselessKeyword, m_props))
molecules = MatchFirst(map(CaselessKeyword, mols))
molecule = MatchFirst(map(CaselessKeyword, mol))
molecule_set = MatchFirst(map(CaselessKeyword, molset))
molecule_sets = MatchFirst(map(CaselessKeyword, molsets))
molecule_identifier = Word(alphas, alphanums + "_" + "[" + "]" + "(" + ")" + "=" + "," + "-" + "+" + "/" + "#" + "@")
INFO_MOLECULES = "\n<soft>To learn more about workspaces, run <cmd>workspace ?</cmd></soft>"
desc = QuotedString("'", escQuote="\\")


def mol_grammar_add(statements, grammar_help):
    """defines the grammar available for managing molecules"""
    statements.append(Forward(add + molecule + molecule_identifier("molecule_identifier"))("add_molecule"))
    grammar_help.append(
        help_dict_create(
            name="add molecule",
            category="Molecules",
            command="add molecule|mol  <name> | <smiles> | <inchi> | <inchkey> | <cid>",
            description="""Adds a given molecule from pubchem to the current working set of molecules. Users can specify a Molecule by Name, a SMILES string, inchi String, Inchkey or its cid.\n
            
            For example:\n 
                - Adding a molecule by name: <cmd> add molecule Aspirin </cmd>\n
                - Adding a molecule by SMILES string: <cmd> add molecule CC(=O)OC1=CC=CC=C1C(=O)O </cmd>\n
                - Adding a molecule by cid: <cmd> add mol 2244 </cmd>\n
                - Adding a molecule by inchikey string: <cmd> add mol  BSYNRYMUTXBXSQ-UHFFFAOYSA-N </cmd>\n
                - Adding a molecule by inchi inchi: <cmd> add mol  InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12) </cmd> \n
                
                Note: when adding a moelcule its name when added to the working list will be the identifying string you added it by. you can use the <cmd> rename </cmd> command to change it""",
            note=INFO_MOLECULES,
        )
    )
    statements.append(Forward(d_isplay + molecule + (molecule_identifier)("molecule_identifier"))("display_molecule"))
    grammar_help.append(
        help_dict_create(
            name="display molecule",
            category="Molecules",
            command="display molecule|mol <name> | <smiles> | <inchi> | <inchkey> |  <cid>",
            description="""Displays a given molecule by first checking the current working set of molecules, then if not in the working set will search for a provided molecule on pubchem. Users can specify a Molecule by Name, a SMILES string, inchi String, Inchkey or its cid.\n
            
            For example:\n 
                - Displaying a molecule by name: <cmd> display molecule Aspirin </cmd>\n
                - Displaying a molecule by SMILES string: <cmd> display molecule CC(=O)OC1=CC=CC=C1C(=O)O </cmd>\n
                - Displaying a molecule by cid: <cmd> display mol 2244 </cmd>\n
                - Displaying a molecule by inchikey string: <cmd> display mol  BSYNRYMUTXBXSQ-UHFFFAOYSA-N </cmd>\n
                - Displaying a molecule by inchi inchi: <cmd> display mol  InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12) </cmd> \n
                """,
            note=INFO_MOLECULES,
        )
    )
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
            command="rename molecule <molecule_identifer_string> name <molecule_name>",
            description="""This command renames a molecule in the current working set to a name you provide it.
            For example:
            I have added a molecule by the molecule 'CC(=O)OC1=CC=CC=C1C(=O)O'to the current working set of molecules, but I want to rename it to 'Aspirin'. the command to do this would be:

            <cmd> rename molecule CC(=O)OC1=CC=CC=C1C(=O)O as Aspirin </cmd>
            """,
            note=INFO_MOLECULES,
        )
    )
    statements.append(
        Forward(
            export
            + molecule
            + (molecule_identifier)("molecule_identifier")
            + Optional(CaselessKeyword("as") + CaselessKeyword("file"))("as_file")
        )("export_molecule")
    )
    grammar_help.append(
        help_dict_create(
            name="export molecule",
            category="Molecules",
            command="export molecule|mol <name> | <smiles> | <inchi> | <inchkey> |  <cid> [as file]",
            description="""exports a molecule from pubchem or the current list to a file named as the molecules given name and or as a dictionary(when in Notebooks) The molecule does not have to be form the current working set, and if not the request will be made to pubchem.
            For Example: 
                    - The following will return a dictionary when called in jupyter notebooks or from the command line it will save it to the current workspace directory as a '.json' file. <cmd> export molecule aspirin <cmd>
                    - The following will  save it to the current workspace directory as a '.json' file. <cmd> export molecule aspirin as file <cmd> """,
            note=INFO_MOLECULES,
        )
    )
    statements.append(Forward(remove + molecule + molecule_identifier("molecule_identifier"))("remove_molecule"))
    grammar_help.append(
        help_dict_create(
            name="remove molecules",
            category="Molecules",
            command="remove molecule|mol <name> | <smiles> | <inchi> | <inchkey> | <formula> | <cid> ",
            description="""removes molecule from the current working set of molecules.
            
            For example:\n 
                - Remove a molecule by name: <cmd> display molecule Aspirin </cmd>\n
                - Remove a molecule by SMILES string: <cmd> display molecule CC(=O)OC1=CC=CC=C1C(=O)O </cmd>\n
                - Remove a molecule by cid: <cmd> display mol 2244 </cmd>\n
                - Remove a molecule by inchikey string: <cmd> display mol  BSYNRYMUTXBXSQ-UHFFFAOYSA-N </cmd>\n
                - Remove a molecule by inchi inchi: <cmd> display mol  InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12) </cmd> \n
                """,
            note=INFO_MOLECULES,
        )
    )
    statements.append(Forward(lister + molecules)("list_molecules"))

    grammar_help.append(
        help_dict_create(
            name="list molecules",
            category="Molecules",
            command="list molecules|mols",
            description="""lists the molecules in the current working set of molecules.
            For example:
            <cmd>list molecules<cmd>
            
            """,
            note=INFO_MOLECULES,
        )
    )

    statements.append(
        Forward(save + molecule_set + a_s + Word(alphas, alphanums + "_")("molecule-set_name"))("save_molecule-set")
    )

    grammar_help.append(
        help_dict_create(
            name="save molecule-set",
            category="Molecules",
            command="save molecule-set|molset as <molecule-set_name>",
            description="""Saves the current molecule working set to the current workspace.
             For example: <cmd> save molecule-set as my_working_set</cmd>""",
            note=INFO_MOLECULES,
        )
    )
    statements.append(
        Forward(load + molecule_set + Word(alphas, alphanums + "_")("molecule-set_name"))("load_molecule-set")
    )

    grammar_help.append(
        help_dict_create(
            name="load molecule-set",
            category="Molecules",
            command="load molecule-set|molset <molecule-set_name>",
            description="""loads the molecules from the current workspace.
            For example: <cmd> load molecule-set my_working_set</cmd>""",
            note=INFO_MOLECULES,
        )
    )

    statements.append(Forward(lister + molecule_sets)("list_molecule-sets"))

    grammar_help.append(
        help_dict_create(
            name="list molecule-sets",
            category="Molecules",
            command="list molecule-sets|molsets",
            description="""lists molecule sets in the current workspace.
            For Example:
            <cmd> List molecule-sets </cmd> or <cmd>list molsets</cmd>""",
            note=INFO_MOLECULES,
        )
    )

    statements.append(Forward(enrich + molecule_set + w_ith + analysis)("load_analysis"))

    grammar_help.append(
        help_dict_create(
            name="enrich molecule-set",
            category="Molecules",
            command="enrich molecule-set with analysis",
            description="Loads the previous results of analysis into the molecule records where the focus of the Analysis was the given molecule",
            note=INFO_MOLECULES,
        )
    )

    statements.append(Forward(clear + analysis + cache)("clear_analysis"))

    grammar_help.append(
        help_dict_create(
            name="clear analysis cache",
            category="Molecules",
            command="clear analysis cache",
            description="Clears the cache of analysis results for the current workspace.",
            note=INFO_MOLECULES,
        )
    )
    statements.append(
        Forward(create + molecule + molecule_identifier("smiles") + name + (Word(alphas, alphanums + "_")("name")))(
            "create_molecule"
        )
    )

    statements.append(Forward(clear + molecules)("clear_molecules"))

    grammar_help.append(
        help_dict_create(
            name="clear Molecules",
            category="Molecules",
            command="clear molecules",
            description="Clears the working set of molecules.",
            note=INFO_MOLECULES,
        )
    )
    statements.append(
        Forward(create + molecule + molecule_identifier("smiles") + name + (Word(alphas, alphanums + "_")("name")))(
            "create_molecule"
        )
    )

    grammar_help.append(
        help_dict_create(
            name="create molecule",
            category="Molecules",
            command="create molecule <smiles_string> name <molecule_name>",
            description="""creates a base molecule and adds it to the current list
            For example: 
                <cmd>create  molecule CC(=O)OC1=CC=CC=C1C(=O)O name my_aspirin</cmd>
                 
                Note it will try and calclculate other idenfiers in the molecule data structure, other propoerties are left as None. """,
            note=INFO_MOLECULES,
        )
    )

    statements.append(
        Forward("@" + molecule_identifier("molecule_identifier") + ">>" + mol_properties("property"))("mol_property")
    )
    description = """allows you to get a molecule property by using one of the propoerties listed below. The molecule can be identified by a name, SMILES, Inchi key or cid.

    For example:
        - Obtain the molecular weight of the molecule known as aspirin. <cmd> @aspirin>>molecular_weight </cmd>
        - Obtain a molecules xlogp value using a SMILES string. <cmd> @CC(=O)OC1=CC=CC=C1C(=O)O>>xlogp </cmd>

     
      Here is a list of valid properties that can be requested. \n""" + str(
        m_props
    )
    grammar_help.append(
        help_dict_create(
            name="@<molecule>",
            category="Molecules",
            command="@(<name> | <smiles> | <inchi> | <inchkey> | <cid>)>><molecule_property_name>",
            description=description,
            note=INFO_MOLECULES,
        )
    )

    statements.append(
        Forward(
            load + molecules + using + file + desc("moles_file") + Optional((merge + w_ith + pubchem))("pubchem_merge")
        )("load_molecules_file")
    )  # From mols file

    grammar_help.append(
        help_dict_create(
            name="load molecules",
            category="Utility",
            command="load molecules using file '<filename>' ",
            description="""load molecules into the molecule working set from a file""",
            note=INFO_MOLECULES,
        )
    )

    statements.append(
        Forward(
            load
            + molecules
            + using
            + CaselessKeyword("dataframe")
            + molecule_identifier("in_dataframe")
            + Optional((merge + w_ith + pubchem))("pubchem_merge")
        )("load_molecules_dataframe")
    )  # From mols file

    grammar_help.append(
        help_dict_create(
            name="load molecules",
            category="Utility",
            command="load molecules using dataframe '<filename>' ",
            description="""load molecules into the molecule working set from a dataframe""",
            note=INFO_MOLECULES,
        )
    )

    statements.append(Forward((export + molecules)("export_molecules")))

    grammar_help.append(
        help_dict_create(
            name="export molecules",
            category="Utility",
            command="export molecules ",
            description="""Exports the molecules in the current Working Set
             
            If the command is issued from a command line the molecule will be exprted to your workspace and named resul_#.csv. # being a incramental number of results sets, with the highest being the latest. 
               
            In Jupyter notebooks the molecules are exported as a pandas Dataframe. """,
            note=INFO_MOLECULES,
        )
    )
