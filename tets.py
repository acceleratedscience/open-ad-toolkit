import pandas
from ad4e_opentoolkit.helpers.output import output_table

lst = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
# df = pandas.DataFrame(lst, columns=['a', 'b', 'c'])
output_table(lst, headers=['a', 'b', 'c'])
