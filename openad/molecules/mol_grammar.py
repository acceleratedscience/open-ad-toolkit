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
) = map(
    CaselessKeyword,
    "get list description using create set unset workspace workspaces context jobs exec\
          as optimize with toolkits toolkit gpu experiment add run save runs show \
              file display history data remove result from inchi inchikey smiles formula name last load results export create rename".split(),
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


def mol_grammar_add(statements, grammar_help):
    """defines the grammar available for managing molecules"""
    statements.append(Forward(add + molecule + molecule_identifier("molecule_identifier"))("add_molecule"))
    grammar_help.append(
        help_dict_create(
            name="add molecule",
            category="Molecules",
            command="add molecule|mol <name> | <smiles> | <inchi> | <inchkey> | <cid>",
            description="adds a given molecule from pubchem to the working List.",
            note=INFO_MOLECULES,
        )
    )
    statements.append(Forward(d_isplay + molecule + (molecule_identifier)("molecule_identifier"))("display_molecule"))
    grammar_help.append(
        help_dict_create(
            name="display molecule",
            category="Molecules",
            command="display molecule|mol <name> | <smiles> | <inchi> | <inchkey> |  <cid>",
            description="Displays a molecule from pubchem or a current working list.",
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
            description="exports a molecule from pubchem or the current list to a file named as the molecules given name and or as a dictionary(when in Notebooks).",
            note=INFO_MOLECULES,
        )
    )
    statements.append(Forward(remove + molecule + molecule_identifier("molecule_identifier"))("remove_molecule"))
    grammar_help.append(
        help_dict_create(
            name="remove molecules",
            category="Molecules",
            command="remove molecule|mol <name> | <smiles> | <inchi> | <inchkey> | <formula> | <cid> ",
            description="removes molecule from working list of molecules.",
            note=INFO_MOLECULES,
        )
    )
    statements.append(Forward(lister + molecules)("list_molecules"))

    grammar_help.append(
        help_dict_create(
            name="list molecules",
            category="Molecules",
            command="list molecules|mols",
            description="lists the molecules in the working list of molecules.",
            note=INFO_MOLECULES,
        )
    )
    statements.append(Forward(add + l_ast + molecule)("add_last_molecule"))
    grammar_help.append(
        help_dict_create(
            name="add last molecule|mol",
            category="Molecules",
            command="add last molecule|mol",
            description="Adds the last molecule retrieved from external source to working list of molecules.",
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
            description="saves a molecule set to the current  workspace.",
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
            description="loads the molecules from the current workspace.",
            note=INFO_MOLECULES,
        )
    )

    statements.append(Forward(lister + molecule_sets)("list_molecule-sets"))

    grammar_help.append(
        help_dict_create(
            name="list molecule-sets",
            category="Molecules",
            command="list molecule-sets|molsets",
            description="lists molecule sets in the current workspace.",
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

    grammar_help.append(
        help_dict_create(
            name="create molecule",
            category="Molecules",
            command="create molecule <smiles_string> name <molecule_name>",
            description="creates a base molecule and adds it to the current list",
            note=INFO_MOLECULES,
        )
    )

    statements.append(
        Forward("@" + molecule_identifier("molecule_identifier") + ">>" + mol_properties("property"))("mol_property")
    )
    description = "allows you to inspect a valid set of properties for a molecule from the following list.\n" + str(
        m_props
    )
    grammar_help.append(
        help_dict_create(
            name="@<molecule>",
            category="Molecules",
            command="@<molecule>>><molecule_property_name>",
            description=description,
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
            description="renames a molecule currently in the working_set",
            note=INFO_MOLECULES,
        )
    )
