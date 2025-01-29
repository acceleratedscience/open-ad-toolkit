"""CONTAINS CONSTANTS AND FUNCTIONS FOR MODELS"""

import platform
import os
from openad.helpers.output import output_error
from langchain_huggingface import HuggingFaceEmbeddings
from langchain_community.embeddings import OllamaEmbeddings
from langchain_community.chat_models import ChatOllama
from genai.schema import TextEmbeddingParameters
from genai import Credentials, Client

# from genai import Model
# from genai.schemas import GenerateParams

from genai import Client, Credentials


# from genai.extensions.langchain.chat_llm import LangChainChatInterface
from genai.extensions.langchain import LangChainEmbeddingsInterface

from genai.extensions.langchain import LangChainInterface


from genai.schema import (
    DecodingMethod,
    ModerationHAP,
    ModerationParameters,
    TextGenerationParameters,
    TextGenerationReturnOptions,
)


# Determine in the sentence transformer embbedings model is installed
# this is currently required for BAM and WATSON
MINI_EMBEDDINGS_MODEL_PRESENT = False

# try:
#    from sentence_transformers import SentenceTransformer
# except:
#    MINI_EMBEDDINGS_MODEL_PRESENT = False


if platform.processor().upper() == "ARM":
    LOCAL_EMBEDDINGS_DEVICE = "mps"
else:
    LOCAL_EMBEDDINGS_DEVICE = "cpu"

LOCAL_MODEL_PATH = "sentence-transformers/all-MiniLM-L6-v2"
LOCAL_MODEL_KWARGS = {"device": LOCAL_EMBEDDINGS_DEVICE}
LOCAL_ENCODE_KWARGS = {"normalize_embeddings": False}

DEFAULT_TELL_ME_MODEL = "OLLAMA"
SUPPORTED_TELL_ME_MODELS = ["BAM", "OLLAMA"]
OLLAMA_HOST = "http://0.0.0.0:11434"

try:
    OLLAMA_HOST = "http://" + os.environ["OLLAMA_HOST"]
except:
    OLLAMA_HOST = "http://0.0.0.0:11434"

SUPPORTED_TELL_ME_MODELS_SETTINGS = {
    "BAM": {
        "model": "ibm/granite-13b-chat-v2",
        "url": "https://bam-api.res.ibm.com",
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


Answer the question based only on the following context: {context}  Question: {question} """,
        "settings": {
            "temperature": 0.3,
            #    "decoding_method": "greedy",
            "max_new_tokens": 2000,
            "min_new_tokens": 1,
            "top_p": 0.3,
            "top_k": 50,
        },
        "embeddings": None,
        "embeddings_api": None,
        "embeddings_model": "sentence-transformers/all-minilm-l6-v2",
    },
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
            # "decoding_method": "greedy",
            "max_new_tokens": 2000,
            "min_new_tokens": 0,
            # "top_p": 0.2,
            # "top_k": 20,
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
        creds = Credentials(api_key=api_key, api_endpoint=SUPPORTED_TELL_ME_MODELS_SETTINGS[service]["url"])

        client = Client(credentials=creds)

        try:
            params = TextGenerationParameters(
                decoding_method=DecodingMethod.GREEDY, max_new_tokens=1536, min_new_tokens=1, temperature=0.3
            )

            model = HuggingFaceEmbeddings(
                client=client,
                model_id=SUPPORTED_TELL_ME_MODELS_SETTINGS[service]["model"],
                parameters=params,
                verbose=True,
            )

            return model, SUPPORTED_TELL_ME_MODELS_SETTINGS[service]["template"]
        except Exception as e:  # pylint: disable=broad-exception-caught
            print(e)
            output_error("Error Loading BAM Model see error Messsage : \n" + e, return_val=False)
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
        ##if MINI_EMBEDDINGS_MODEL_PRESENT is False:
        ##    return False
        try:
            creds = Credentials(
                api_key=api_key,
                api_endpoint=SUPPORTED_TELL_ME_MODELS_SETTINGS[service]["url"],
            )
            client = Client(credentials=creds)
            embeddings = LangChainEmbeddingsInterface(
                client=client,
                model_id=SUPPORTED_TELL_ME_MODELS_SETTINGS[service]["embeddings_model"],
                parameters=TextEmbeddingParameters(truncate_input_tokens=True),
            )

        except Exception as e:
            print(e)
            print("Error: cannot initialise embeddings, check BAM requirements isntalled")
            raise Exception("Error: cannot initialise embeddings, check BAM requirements isntalled") from e

    return embeddings
