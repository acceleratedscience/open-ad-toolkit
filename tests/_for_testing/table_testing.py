"""
This files contains a number of different table types meant to test output_table().
"""

import os
import sys
import pandas as pd

parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, parent_dir)
from openad.helpers.output import output_table


# fmt: off
headers_1 = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T"]
headers_2 = ["A\nA\nA", "B\nB\nB", "C\nC\nC", "D\nD\nD", "E\nE\nE", "F\nF\nF", "G\nG\nG", "H\nH\nH", "I\nI\nI", "J\nJ\nJ", "K\nK\nK", "L\nL\nL", "M\nM\nM", "N\nN\nN", "O\nO\nO", "P\nP\nP", "Q\nQ\nQ", "R\nR\nR", "S\nS\nS", "T\nT\nT"]
# fmt: on

# List, 100 rows
table_list_long = []
for i in range(1, 101):
    table_list_long.append(("%03d" % i,) * 20)

# List, 20 rows, multi-line line header
table_list_short = []
for i in range(1, 11):
    table_list_short.append(("%02d" % i,) * 20)

# Dataframe, 20 rows
table_df_short = pd.DataFrame(table_list_short)

# Dataframe, 100 rows
table_df_long = pd.DataFrame(table_list_long)

# Dataframe, 20 rows, embedded columns
table_df_short_columns = pd.DataFrame(
    table_list_short,
    columns=["D", "A", "T", "A", "F", "R", "A", "M", "E", "-", "C", "O", "L", "U", "M", "N", "S", "-", "O", "K"],
)

# Dataframe, 100 rows, embedded columns
table_df_long_columns = pd.DataFrame(
    table_list_long,
    columns=["D", "A", "T", "A", "F", "R", "A", "M", "E", "-", "C", "O", "L", "U", "M", "N", "S", "-", "O", "K"],
)


print("++++++++++")

#
#
# SHORT LIST TABLE

# Short list-table
# output_table(table_list_short, headers=headers_1, _display="terminal")

# Short list-table, multi-line headers
# output_table(table_list_short, headers=headers_2, _display="terminal")

# Short list-table, no follow-up commands
# output_table(table_list_short, headers=headers_1, is_data=False, _display="terminal")

# Short list-table, custom padding
# output_table(table_list_short, headers=headers_1, pad=3, _display="terminal")

# Short list-table, custom top padding
# output_table(table_list_short, headers=headers_1, pad_top=3, _display="terminal")

# Short list-table, custom bottom padding
# output_table(table_list_short, headers=headers_1, pad_btm=3, _display="terminal")

#
#
# PAGINATED LIST TABLE

# Paginated list-table
# output_table(table_list_long, headers=headers_1, _display="terminal")

# Paginated list-table, multi-line header
# output_table(table_list_long, headers=headers_2, _display="terminal")

# Paginated list-table w/o follow-up commands
# output_table(table_list_long, headers=headers_1, is_data=False, _display="terminal")

#
#
# SHORT DATAFRAME TABLE

# Short dataframe-table
# output_table(table_df_short, headers=headers_1, _display="terminal")

# Short dataframe-table, no headers
# output_table(table_df_short, _display="terminal")

# Short dataframe-table, headers from dataframe
# output_table(table_df_short_columns, _display="terminal")

# Short dataframe-table, styled
# (Has no effect in the terminal, but is supported for Jupyter)
output_table(table_df_short.style.set_properties(**{"text-align": "left"}), headers=headers_1, _display="terminal")

#
#
# PAGINATED DATAFRAME TABLE

# Paginated dataframe-table
# output_table(table_df_long, headers=headers_2, _display="terminal")

# Paginated dataframe-table, no headers
# output_table(table_df_long, _display="terminal")

# Paginated dataframe-table, headers from dataframe
# output_table(table_df_long_columns, _display="terminal")

print("++++++++++")
