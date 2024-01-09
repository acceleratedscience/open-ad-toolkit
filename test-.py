from tabulate import tabulate
from openad.helpers.output import _paginated_output

headers = ["X\nX", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R"]
data = [(1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000) * 2]
data = data * 5
table = tabulate(data, headers=headers, tablefmt="simple", showindex=False, numalign="left")

_paginated_output(table, headers=headers)

print(table)
