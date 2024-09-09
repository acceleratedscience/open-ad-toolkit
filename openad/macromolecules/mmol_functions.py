import copy
import requests
from Bio import Entrez
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# NCBI identification
# https://biopython.org/docs/latest/Tutorial/chapter_entrez.html#entrez-guidelines
Entrez.tool = "IBM Research - OpenAD"
Entrez.email = "phil.downey1@ibm.com"  # Email required by NCBI
# Entrez.api_key = ""  # Allows 10 queries/s instead of 3 queries/s - See https://tinyurl.com/ncbi-api-key


# Eg. SATVSEINSETDFVAKNDQFIALTKDTTAHIQSNSLQSVEELHSSTINGVKFEEYLKSQIATIGENLVVRRFATLKAGANGVVNGYIHTNGRVGVVIAAACDSAEVASKSRDLLRQICMH
def get_protein(identifier):
    "Retrieve protein by its identifier"

    # For now FASTA is the only identifier type.
    fasta_string = identifier

    try:
        # Create a Seq object
        protein_seq = Seq(fasta_string)

        # Analyze the protein sequence
        protein_analysis = ProteinAnalysis(str(protein_seq))

        # Get protein properties
        properties = {
            "molecular_weight": protein_analysis.molecular_weight(),
            "aromaticity": protein_analysis.aromaticity(),
            "instability_index": protein_analysis.instability_index(),
            "isoelectric_point": protein_analysis.isoelectric_point(),
            "secondary_structure_fraction": protein_analysis.secondary_structure_fraction(),
        }

        # Store biopython as the source of the properties
        property_sources = {}
        for key in properties:
            property_sources[key] = "biopython"

        # Compile result
        result = {
            "identifiers": {
                "fasta": fasta_string,
            },
            "synonyms": [],
            "properties": properties,
            "property_sources": property_sources,
        }

        return result

    except Exception as e:
        print(e)
        return False


# Eg. 4PC7
def get_protein2(identifier):
    print("get_protein", identifier)

    # Search for the protein by identifier
    search_handle = Entrez.esearch(db="protein", term=identifier)
    search_results = Entrez.read(search_handle)
    search_handle.close()

    if not search_results["IdList"]:
        print("No matching protein entries found.")
        return

    # Get the first matching protein ID
    protein_id = search_results["IdList"][0]

    print("protein_id", protein_id)

    # Fetch the protein data
    fetch_handle = Entrez.efetch(db="protein", id=protein_id, rettype="gb", retmode="text")
    print(1)
    protein_data = fetch_handle.read()
    print(2)
    fetch_handle.close()
    print(3)

    print("protein_data", protein_data)

    return protein_data


def get_protein3(fasta_sequence):
    print(1)

    # Perform BLAST search
    result_handle = NCBIWWW.qblast("blastp", "pdb", fasta_sequence)

    print(2)

    # Parse BLAST results
    blast_record = NCBIXML.read(result_handle)

    print(3)

    if not blast_record.alignments:
        print("No matching PDB entries found.")
        return

    # Get the first matching PDB ID
    top_hit = blast_record.alignments[0]
    pdb_id = top_hit.accession

    print(4)

    # Fetch the PDB data
    fetch_handle = Entrez.efetch(db="structure", id=pdb_id, rettype="pdb", retmode="text")
    print(5)
    pdb_data = fetch_handle.read()
    print(6)
    fetch_handle.close()
    print(7)

    print("pdb_data", pdb_data)

    return pdb_data


# Not great using web api
def fetch_pdb_data(fasta_string):

    # Search for the sequence in the PDB
    search_url = f'https://search.rcsb.org/rcsbsearch/v2/query?json={{"query":{{"type":"terminal","service":"sequence","parameters":{{"value":"{fasta_string}","evalue_cutoff":0.1,"identity_cutoff":0.95}}}},"return_type":"entry"}}'
    search_response = requests.get(search_url)
    search_results = search_response.json()

    if not search_results.get("result_set"):
        return "No matching PDB entries found."

    # Get the first matching PDB ID
    pdb_id = search_results["result_set"][0]["identifier"]

    # Fetch the PDB data
    pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    pdb_response = requests.get(pdb_url)

    if pdb_response.status_code != 200:
        return "Failed to retrieve PDB data."

    pdb_data = pdb_response.text

    return pdb_data
