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
)

(protein, show) = map(
    CaselessKeyword,
    "protein show".split(),
)
protein_identifier = Word(alphas + "-" + "*")
desc = QuotedString("'", escQuote="\\")


def prot_grammar_add(statements, grammar_help):
    """
    Defines the grammar available for managing proteins.
    """

    # ---
    # Show individual protein in browser.
    statements.append(
        Forward(show("show") + protein + (protein_identifier | desc)("protein_identifier"))("show_protein")
    )
    grammar_help.append(
        help_dict_create(
            name="show protein",
            category="Proteins",
            command="show protein <fasta>",
            description="""
Inspect a protein in the browser.

Examples:
- <cmd>show protein MAAVLLFLLVPGAGLAMRLLGLLLVGLPV</cmd>
""",
        )
    )
