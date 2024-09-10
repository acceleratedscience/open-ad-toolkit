import json
import requests
from Bio import Entrez
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from openad.helpers.general import encode_uri_component


def mmol_from_identifier(identifier):
    """
    Fetch a macromolecule from an identifier.

    Parameters:
        identifier: str
            FASTA or PDB id.

    Returns:
        success: bool
            Whether the request was successful.
        mmol_data: dict
            The macromolecule data.
    """

    # Try fetching the mmol by its PDB ID.
    if len(identifier) == 4:
        success, cif_data = fetch_pdb_file(identifier)

    # Try fetching the mmol by a FASTA sequence.
    if not success:
        success, cif_data = search_fasta_sequence(identifier)

    # Fail
    if not success:
        return False, "Failed to retrieve macromolecule data."

    return True, cif_data


def search_fasta_sequence(fasta_string, sequence_type="protein", return_first=True):
    """
    Search the RCSB PDB for a given FASTA sequence.

    Parameters:
        fasta_string: str
            The FASTA sequence to search for.
        sequence_type: 'protein' | 'dna' | 'rna'
            The type of sequence to search for.
            https://search.rcsb.org/#search-services

    Returns:

    """

    # https://search.rcsb.org/#search-example-3
    search_dict = {
        "query": {
            "type": "terminal",
            "service": "sequence",
            "parameters": {
                "evalue_cutoff": 1,
                "identity_cutoff": 0.9,
                "sequence_type": sequence_type,
                "value": fasta_string,
            },
        },
        "request_options": {"scoring_strategy": "sequence"},
        "return_type": "entry",
    }

    # Make search_dict URL safe
    search_str = json.dumps(search_dict)
    search_str = search_str.replace(" ", "")
    search_str = encode_uri_component(search_str)

    # Search for the sequence in the PDB
    search_url = f"https://search.rcsb.org/rcsbsearch/v2/query?json={search_str}"
    search_response = requests.get(search_url)
    search_results = search_response.json() if search_response.status_code == 200 else {}

    # Error handling
    if not search_results.get("result_set"):
        return False, "No matching PDB entries found."

    # Return the first result if it's a 100% match
    if len(search_results["result_set"]) > 0 and search_results["result_set"][0].get("score") == 1:
        pdb_id = search_results["result_set"][0]["identifier"]
        return fetch_pdb_file(pdb_id)
    else:
        return False, search_results["result_set"]


def fetch_pdb_file(pdb_id, file_format="cif"):
    """
    Fetch a PDB file by its ID from rscb.org.

    Parameters:
        pdb_id: str
            The PDB ID to fetch.
        file_format: 'cif' | 'pdb'
            The format of the PDB file to fetch.

    Returns:
        success: bool
            Whether the request was successful.
        file_data: str
            The content of the (cif) file.
    """

    # Fetch the file data.
    pdb_url = f"https://files.rcsb.org/download/{pdb_id}.{file_format}"
    pdb_response = requests.get(pdb_url)

    if pdb_response.status_code != 200:
        return False, "Failed to retrieve PDB data."

    file_data = pdb_response.text
    return True, file_data


def ncbi_search(identifier):
    """
    Search the NCBI database.

    Test with:
    ncbi_search("P0A9Q1")

    Currently not used.
    """

    # NCBI identification
    # https://biopython.org/docs/latest/Tutorial/chapter_entrez.html#entrez-guidelines
    Entrez.tool = "IBM Research - OpenAD"
    Entrez.email = "phil.downey1@ibm.com"  # Email required by NCBI
    # Entrez.api_key = ""  # Allows 10 queries/s instead of 3 queries/s - See https://tinyurl.com/ncbi-api-key

    # Search for the protein by identifier
    search_handle = Entrez.esearch(db="protein", term=identifier)
    search_results = Entrez.read(search_handle)
    search_handle.close()

    if not search_results["IdList"]:
        print("No matching protein entries found.")
        return

    # Get the first matching protein ID
    protein_id = search_results["IdList"][0]

    # Fetch the protein data
    fetch_handle = Entrez.efetch(db="protein", id=protein_id, rettype="gb", retmode="text")
    protein_data = fetch_handle.read()
    fetch_handle.close()

    print("protein_data", protein_data)

    return protein_data


def parse_cif_block(cif_block):
    """
    Parse a gemmi CIF block into a dictionary.

    Parameters:
        cif_block: a gemmi CIF block.

    Returns:
        data: dict
            The parsed CIF dictionary
    """

    data = {}
    for item in cif_block:
        if item.pair is not None:
            cat = item.pair[0].split(".")[0].lstrip("_")
            key = item.pair[0].split(".")[1]
            if cat not in data:
                data[cat] = {}
            data[cat][key] = item.pair[1]

    return data


# fmt: off
# For testing
if __name__ == "__main__":
    # x, y = search_fasta_sequence("IINVKTSLKTIIKNALDKIQX") # Positive result
    # x, y = search_fasta_sequence("MSKGEELFTTYQDKDTAGHKHYGSHQYAERVGGMPEYMFTQVTGDRCDNAQYNGVLYQWDAMKKYGGERQGIVQLKPGTFGAVK") # No results
    # print(x, y)
    # ncbi_search("P0A9Q1")
    pass
