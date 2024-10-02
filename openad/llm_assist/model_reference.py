"""CONTAINS CONSTANTS AND FUNCTIONS FOR MODELS"""

import platform
import os
from openad.helpers.output import output_error
from langchain_openai import ChatOpenAI

from langchain_openai import OpenAIEmbeddings
from langchain_community.embeddings import HuggingFaceEmbeddings
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
# SUPPORTED_TELL_ME_MODELS = ["BAM", "OPENAI", "OLLAMA"]
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
            "temperature": 0.5,
            "decoding_method": "greedy",
            "max_new_tokens": 5000,
            "min_new_tokens": 1,
            "top_p": 0.3,
            "top_k": 50,
        },
        "embeddings": None,
        "embeddings_api": None,
    },
    "OLLAMA": {
        "model": "llama3.1:8b",
        # "model": "instructlab/granite-7b-lab:latest",
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


Answer the question based only on the following context: {context}  Question: {question} """,
        "template1": """You are a technical documentation writer and when responding follow the following rules:
                - Respond like you were writing a refernce guide for a software package
                - Format All Command Syntax, Clauses, Examples or Option  Syntax in codeblock ipython Markdown
                - Format all Command Syntax, Options or clause quotations in codeblock ipython Markdown
                - Only format codeblocks one line at a time and place them  on single lines
                - For each instruction used in an answer also provide full command syntax with clauses and options in codeblock format. for example " Use the `search collection` with the 'PubChem' collection to search for papers and molecules.   \n\n command: ` search collection '<collection name or key>' for '<search string>' using ( [ page_size=<int> system_id=<system_id> edit_distance=<integer> display_first=<integer>]) show (data|docs) [ estimate only|return as data|save as '<csv_filename>' ] ` \n
                \n For Example: ` search collection 'PubChem' for 'Ibuprofen' show ( data ) ` \n"
                - Provide All syntax, clauses, Options, Parameters and Examples separated by "\n" for a command when answering a question with no leading spaces on the line
                - Compounds and Molecules are the same concept
                - Always explain using the full name not short form of a name
                - Never refer to source files from the embeddings
                - after explaning a command  tell them how to go to the help using `<command> ?` substituing the command into the string
                - if a Property is mentioned look for commands that mention the Property exactly
                - respond with a output format as per following example
                - Commands in the documentation are bracketed in the <cmd>  </cmd> tags
                - prioritise command responses from examples in <cmd> tags
                '''Command: <put command syntax here >

                Description: <brief description of funciton

                Parameters: < Tell the user what Parameters are available for the comand>

                Examples:  < examples of how to use the function> '''

Answer the question based only on the following 
context: {context} 

Question: {question}  

Answer:""",
        "settings": {
            "temperature": 0.2,
            "decoding_method": "greedy",
            "max_new_tokens": 5000,
            "min_new_tokens": 1,
            "top_p": 0.3,
            "top_k": 50,
        },
        "embeddings": None,
        "embeddings_api": None,
    },
    "OPENAI": {
        "model": "gpt-3.5-turbo",
        "url": None,
        "template": """When responding follow the following rules:
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
        
            Answer the question based only on the following 
            
            context: {context} 
             
            Question: {question}.  Please provide information on options in response. Please format all code syntax, clauses and options in codeblock formatting on single lines in the response. 
            
            Answer:""",
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
            output_error("Error Loading  Model see error Messsage : \n" + e, return_val=False)
            return None, None
    elif service == "OLLAMA":
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
    elif service == "OLLAMA":
        try:
            embeddings = OllamaEmbeddings(
                model="pankajrajdeo/sentence-transformers_all-minilm-l6-v2",
                # model="jmorgan/all-minilm-l6",
                base_url=OLLAMA_HOST,
                # model_kwargs={"truncation": True},
            )
            # embeddings = OllamaEmbeddings(model="all-minilm", base_url=OLLAMA_HOST)
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
                model_id="sentence-transformers/all-minilm-l6-v2",
                parameters=TextEmbeddingParameters(truncate_input_tokens=True),
            )

        # embeddings = HuggingFaceEmbeddings(
        #    model_name=LOCAL_MODEL_PATH,
        #    model_kwargs=LOCAL_MODEL_KWARGS,
        #    encode_kwargs=LOCAL_ENCODE_KWARGS,
        # )

        except Exception as e:
            print(e)
            print("Error: cannot initialise embeddings, check BAM requirements isntalled")
            raise Exception("Error: cannot initialise embeddings, check BAM requirements isntalled") from e

    """elif service == "WATSONX":
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
        # If not refreshing the database, check to see if the database exists"""
    return embeddings
