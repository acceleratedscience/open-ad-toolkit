"""CONTAINS CONSTANTS AND FUNCTIONS FOR MODELS"""
import platform

from openad.helpers.output import output_error
from langchain.chat_models import ChatOpenAI
from langchain.embeddings.openai import OpenAIEmbeddings
from langchain.embeddings import HuggingFaceEmbeddings

from genai import Credentials, Client

# from genai import Model
# from genai.schemas import GenerateParams
from genai import Client, Credentials
from genai.extensions.langchain.chat_llm import LangChainChatInterface
from genai.extensions.langchain import LangChainInterface
from genai.text.generation import (
    DecodingMethod,
    ModerationHAP,
    ModerationParameters,
    TextGenerationParameters,
    TextGenerationReturnOptions,
)

# Determine in the sentence transformer embbedings model is installed
# this is currently required for BAM and WATSON
MINI_EMBEDDINGS_MODEL_PRESENT = True
try:
    from sentence_transformers import SentenceTransformer
except:
    MINI_EMBEDDINGS_MODEL_PRESENT = False


if platform.processor().upper() == "ARM":
    LOCAL_EMBEDDINGS_DEVICE = "mps"
else:
    LOCAL_EMBEDDINGS_DEVICE = "cpu"

LOCAL_MODEL_PATH = "sentence-transformers/all-MiniLM-L6-v2"
LOCAL_MODEL_KWARGS = {"device": LOCAL_EMBEDDINGS_DEVICE}
LOCAL_ENCODE_KWARGS = {"normalize_embeddings": False}

DEFAULT_TELL_ME_MODEL = "OPENAI"
SUPPORTED_TELL_ME_MODELS = ["BAM", "OPENAI"]

SUPPORTED_TELL_ME_MODELS_SETTINGS = {
    "BAM": {
        "model": "ibm/granite-20b-code-instruct-v1",
        "model2": "mistralai/mistral-7b-instruct-v0-2",
        "url": "https://bam-api.res.ibm.com",
        "template": """  When responding follow the following rules:
                - Answer and format like a Technical Documentation writer concisely and to the point
                - Format All Command Syntax, Clauses, Examples or Option Syntax in codeblock Markdown
                - Format all Command Syntax, Options or clause quotations in codeblock Markdown
                - Only format codeblocks one line at a time and place them  on single lines
                - For each instruction used in an answer also provide full command syntax with clauses and options in codeblock format. for example " Use the `search collection` with the 'PubChem' collection to search for papers and molecules.   \n\n command: ` search collection '<collection name or key>' for '<search string>' using ( [ page_size=<int> system_id=<system_id> edit_distance=<integer> display_first=<integer>]) show (data|docs) [ estimate only|return as data|save as '<csv_filename>' ] ` \n
                \n For Example: ` search collection 'PubChem' for 'Ibuprofen' show ( data ) ` \n"
                - Provide All syntax, clauses, Options, Parameters and Examples separated by "\n" for a command when answering a question with no leading spaces on the line
                - ensure bullet lines are indented consistently
                - Compounds and Molecules are the same concept
                - smiles or inchi strings are definitions of compounds or smiles
                - Always explain using the full name not short form of a name
                - do not repeat instructions unless necessary in answers
                - do not return data in a table format


       

Answer the question based only on the following context: {context}  Question: {question} """,
        "settings": {"temperature": 0.5, "decoding_method": "greedy", "max_new_tokens": 1536, "min_new_tokens": 1},
        "embeddings": None,
        "embeddings_api": None,
    },
    "OPENAI": {
        "model": "gpt-3.5-turbo",
        "url": None,
        "template": """  When responding follow the following rules:
            - Respond like a technical helpful writer
            - Respond succinctly wihout repetition
            - Explain what any requested commands do
            - Provide All syntax, Options, Parameters and Examples when answering a question
            - Provide All Command Syntax, Clauses or Option Syntax in codeblock Markdown
            - Only format one line at a time. codeblocks should not go over lines
            - No "\n" characters in codeblocks
            - Use this correct version of an example command in codeblock format ``` search collection 'PubChem' for 'Ibuprofen' show (data)  ``` 
            - Always format the answer 
            - do not return data in a table format
        
            Answer the question based only on the following context: {context} 
             
            Question: {question}.  Please provide information on options in response. Please format all code syntax, clauses and options in codeblock formatting on single lines in the response.""",
        "settings": {"temperature": None, "decoding_method": "greedy", "max_new_tokens": 1536, "min_new_tokens": 1},
        "embeddings": "text-embedding-ada-002",
        "embeddings_api": None,
    },
}


def get_tell_me_model(service: str, api_key: str):
    """gets the Model Object and Template for the given Tell Me Model Service"""
    if service not in SUPPORTED_TELL_ME_MODELS:
        return None, None

    if service == "OPENAI":
        try:
            model = ChatOpenAI(
                model_name=SUPPORTED_TELL_ME_MODELS_SETTINGS[service]["model"],
                openai_api_key=api_key,
            )
            return model, SUPPORTED_TELL_ME_MODELS_SETTINGS[service]["template"]
        except Exception as e:  # pylint: disable=broad-exception-caught
            output_error("Error Loading OPENAI Model see error Messsage : \n" + e, return_val=False)
            return None, None

    elif service == "BAM":
        if MINI_EMBEDDINGS_MODEL_PRESENT is False:
            output_error(
                "Error: Loading BAM Model you need to install `sentence-transformers` : \n" + e, return_val=False
            )
            return False

        creds = Credentials(api_key=api_key, api_endpoint=SUPPORTED_TELL_ME_MODELS_SETTINGS[service]["url"])

        client = Client(credentials=creds)

        try:
            params = TextGenerationParameters(
                decoding_method=DecodingMethod.GREEDY, max_new_tokens=1536, min_new_tokens=1, temperature=0.3
            )
            # model = Model(
            # model=SUPPORTED_TELL_ME_MODELS_SETTINGS[service]["template"], credentials=creds, params=params
            # )
            model = LangChainInterface(
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
    if service == "OPENAI":
        try:
            embeddings = OpenAIEmbeddings(
                openai_api_key=api_key,
            )
        except Exception as e:
            raise Exception(
                "Error: cannot initialise embeddings, check API Key"
            ) from e  # pylint: disable=broad-exception-raised
    elif service == "BAM":
        if MINI_EMBEDDINGS_MODEL_PRESENT is False:
            return False
        try:
            # creds = Credentials(
            #    api_key=api_key,
            #    api_endpoint=SUPPORTED_TELL_ME_MODELS_SETTINGS[service]["url"],
            #    model_id=LOCAL_MODEL_PATH,
            # )
            # embeddings = Client(credentials=creds)

            embeddings = HuggingFaceEmbeddings(
                model_name=LOCAL_MODEL_PATH,
                model_kwargs=LOCAL_MODEL_KWARGS,
                encode_kwargs=LOCAL_ENCODE_KWARGS,
            )

        except Exception as e:
            print(e)
            print("Error: cannot initialise embeddings, check BAM requirements isntalled")
            raise Exception("Error: cannot initialise embeddings, check BAM requirements isntalled") from e

    elif service == "WATSONX":
        if MINI_EMBEDDINGS_MODEL_PRESENT is False:
            return False
        try:
            embeddings = HuggingFaceEmbeddings(
                model_name=LOCAL_MODEL_PATH,
                model_kwargs=LOCAL_MODEL_KWARGS,
                encode_kwargs=LOCAL_ENCODE_KWARGS,
            )

        except Exception as e:
            raise Exception(
                "Error: cannot initialise embeddings, check API Key"
            ) from e  # pylint: disable=broad-exception-raised
        # If not refreshing the database, check to see if the database exists
    return embeddings
