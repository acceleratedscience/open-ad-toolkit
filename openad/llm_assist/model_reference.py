"""CONTAINS CONSTANTS AND FUNCTIONS FOR MODELS"""

import platform
import os
from openad.helpers.output import output_error


# from langchain_huggingface import HuggingFaceEmbeddings
from langchain_community.embeddings import OllamaEmbeddings
from langchain_community.chat_models import ChatOllama


# Determine in the sentence transformer embbedings model is installed
# this is currently required for BAM and WATSON
MINI_EMBEDDINGS_MODEL_PRESENT = False


if platform.processor().upper() == "ARM":
    LOCAL_EMBEDDINGS_DEVICE = "mps"
else:
    LOCAL_EMBEDDINGS_DEVICE = "cpu"

# Local Not currently used but leaving open for local docling
LOCAL_MODEL_PATH = "sentence-transformers/all-MiniLM-L6-v2"
LOCAL_MODEL_KWARGS = {"device": LOCAL_EMBEDDINGS_DEVICE}
LOCAL_ENCODE_KWARGS = {"normalize_embeddings": False}

DEFAULT_TELL_ME_MODEL = "OLLAMA"
SUPPORTED_TELL_ME_MODELS = ["OLLAMA"]
OLLAMA_HOST = "http://0.0.0.0:11434"

try:
    OLLAMA_HOST = "http://" + os.environ["OLLAMA_HOST"]
except:
    OLLAMA_HOST = "http://0.0.0.0:11434"

SUPPORTED_TELL_ME_MODELS_SETTINGS = {
    "OLLAMA": {
        "model": "granite3.1-dense:8b-instruct-q4_1",
        "url": OLLAMA_HOST,
        "template": """  When responding follow the following rules:
                - Answer and format like a Technical Documentation writer concisely and to the point
                - Format All Command Syntax, Clauses, Examples or Option Syntax in codeblock ipython Markdown
                - Format all Command Syntax, Options or clause quotations in codeblock ipython Markdown
                - Only format codeblocks one line at a time and place them  on single lines
                - For each instruction used in an answer also provide full command syntax with clauses and options in codeblock format. for example " Use the `search collection` with the 'PubChem' collection to search for papers and molecules.   \n\n command: ` search collection '<collection name or key>' for '<search string>' using ( [ page_size=<int> system_id=<system_id> edit_distance=<integer> display_first=<integer>]) show (data|docs) [ estimate only|return as data|save as '<csv_filename>' ] ` \n
                \n For Example: ` search collection 'PubChem' for 'Ibuprofen' show ( data ) ` \n"
                - Provide All syntax, clauses, Options, Parameters and Examples separated by "\n" for a command when answering a question with no leading spaces on the line
                - ensure bullet lines are indented consistently
                - Compounds and Molecules are the same concept
                - smiles or inchi strings are definitions of compounds or smiles
                - Always explain using the full name not short form of a name
                - do not return data in a table format
                - Always list all options and parameters that are documented for a command
                - All ways keep syntax and options provided for a command to be sourced from a single document do not mix
                - do not provide examples of a command not shown in examples for a command, and use exact syntax (include clauses and brackets)
                - if the question asks for "list all paramaters" or "list all options" for a command display both the Optional Parameters and Required Parameters for the requested command
                
Answer the question based only on the following context using provided embeddings only: {context}  Question: {question} """,
        "settings": {
            "temperature": 0.4,
            "max_new_tokens": 2000,
            "min_new_tokens": 0,
            "presence_penalty": 0,
            "frequency_penalty": 0,
        },
        "embeddings": None,
        "embeddings_api": None,
        "embeddings_model": "all-minilm:33m",
    },
}


def get_tell_me_model(service: str, api_key: str):
    """gets the Model Object and Template for the given Tell Me Model Service"""
    # api_key used if needed for authentication
    if service not in SUPPORTED_TELL_ME_MODELS:
        return None, None

    if service == "OLLAMA":
        try:
            model = ChatOllama(model=SUPPORTED_TELL_ME_MODELS_SETTINGS[service]["model"], base_url=OLLAMA_HOST)

            return model, SUPPORTED_TELL_ME_MODELS_SETTINGS[service]["template"]
        except Exception as e:  # pylint: disable=broad-exception-caught
            output_error("Error Loading  Model see error Messsage : \n" + e, return_val=False)
            return None, None
    elif service == "BAM":
        # This is No Longer supported

        return None, None


def get_embeddings_model(service: str, api_key: str):
    """get the embeddings model based on the selected Supported Service"""
    embeddings = None
    if service == "OLLAMA":
        try:
            embeddings = OllamaEmbeddings(
                model=SUPPORTED_TELL_ME_MODELS_SETTINGS[service]["embeddings_model"],
                base_url=OLLAMA_HOST,
                model_kwargs={"truncation": True},
            )
        except Exception as e:
            raise Exception(
                "Error: cannot initialise embeddings, check API Key"
            ) from e  # pylint: disable=broad-exception-raised
    elif service == "BAM":
        # BAM Not supported
        return None

    return embeddings
