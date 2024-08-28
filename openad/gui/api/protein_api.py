"""
Protein API
"""

import os
import json
import shutil
from urllib.parse import unquote
from rdkit import Chem
from flask import Response, request

from openad.macromolecules.mmol_functions import get_protein


from openad.helpers.files import open_file
from openad.helpers.output import output_error, output_table
from openad.helpers.json_decimal_encoder import DecimalEncoder


class ProteinApi:
    """
    All the API endpoints related to proteins.
    The API endpoints are called from gui_routes.py.
    """

    def __init__(self, cmd_pointer):
        self.cmd_pointer = cmd_pointer

    def get_prot_data(self):
        """
        Get protein data.
        Used when requesting a protein by its FASTA identifier.
        """

        data = json.loads(request.data) if request.data else {}
        identifier = data["identifier"] if "identifier" in data else ""

        if not identifier:
            response = Response(None, status=500)
            response.status = "No identifier provided."
            return response

        prot = get_protein(identifier)

        # Fail
        if not prot:
            response = Response(None, status=500)
            response.status = f"No protein found with provided identifier '{identifier}'"
            return response

        # Success
        else:
            return prot, 200
