import pyparsing as py

#############################################
# region - General

# Optional quoted string
quoted = py.QuotedString("'", escQuote="\\")
# quoted.setParseAction(py.removeQuotes)
opt_quoted = py.Optional(quoted)


# endregion
#############################################
# region - Lists

# Open and close square brackets, with optional space on inside
sb_open = py.Suppress(py.Literal("[") + py.Optional(py.White()))
sb_close = py.Suppress(py.Optional(py.White()) + py.Literal("]"))

# Comma with optional spaces before/after
comma = py.Suppress(py.Optional(py.White()) + py.Literal(",") + py.Optional(py.White()))


# endregion
#############################################
# region - Molecules

# Molecule(s) keywords
molecule = py.CaselessKeyword("molecule") | py.CaselessKeyword("mol")
molecules = py.CaselessKeyword("molecules") | py.CaselessKeyword("mols")

# Recursive molecule identifier that allows to parse a square brackets list
# with molecule identifiers that contain square brackets themselves.
# This requires that the list is the end of the command.
molecule_identifier = py.Forward()
_mol_chars = py.alphanums + "_()=-+/\\#@.,*;"
_mol_str = py.Word(_mol_chars)
_mol_str_brackets = py.Combine(py.Literal("[") + molecule_identifier + py.Literal("]"))
molecule_identifier <<= py.Combine(py.OneOrMore(_mol_str | _mol_str_brackets))

# List of identifier strings separated by comas
_delimited_molecule_identifiers = py.delimitedList(molecule_identifier, delim=comma)

# List of strings separated by comas and encapsulated in square brackets
molecule_identifier_list = sb_open + _delimited_molecule_identifiers + sb_close


# endregion
