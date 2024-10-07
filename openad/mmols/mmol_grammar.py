from openad.core.help import help_dict_create
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
)

(add, show) = map(
    CaselessKeyword,
    "add show".split(),
)

_mmol = ["macromolecule", "macromol", "mmol", "protein", "prot"]
mmol = MatchFirst(map(CaselessKeyword, _mmol))
mmol_identifier = Word(alphas + "-" + "*")
desc = QuotedString("'", escQuote="\\")


def mmol_grammar_add(statements, grammar_help):
    """
    Defines the grammar available for managing macromolecules.
    """

    # ---
    # Show individual protein in browser.
    statements.append(Forward(show("show") + mmol + (mmol_identifier | desc)("mmol_identifier"))("show_mmol"))
    grammar_help.append(
        help_dict_create(
            name="show macromolecule",
            category="Macromolecules",
            command="show mmol|protein <fasta> | '<pdb_id>'",
            description=f"""Launch the molecule viewer { 'in your browser ' if is_notebook_mode() else '' }to visualize your macromolecule and inspect its properties.

Examples:
- Show a protein by its PDBe ID:
  <cmd>show mmol '2g64'</cmd>

- Show a protein by its FASTA string:
  <cmd>show protein MAKWVCKICGYIYDEDAGDPDNGISPGTKFEELPDDWVCPICGAPKSEFEKLED</cmd>
""",
        )
    )
