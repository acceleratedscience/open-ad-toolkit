"""
Functions called by the macromolecule commands.
"""

import glob
import pickle
import os
import shutil
import re
import json
import urllib.parse
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

# Helpers
from openad.helpers.general import confirm_prompt
from openad.helpers.output import output_text, output_table, output_warning, output_error, output_success
from openad.helpers.output_msgs import msg
from openad.helpers.format_columns import single_value_columns, name_and_value_columns

# Protein functions
# from openad.macromolecules.mmol_functions import show_protein


def show_protein(cmd_pointer, inp):
    from openad.gui.gui_launcher import gui_init

    protein_identifier = inp.as_dict()["protein_identifier"]
    path = "prot/" + urllib.parse.quote(protein_identifier, safe="")
    gui_init(cmd_pointer, path)
