""" this library Automates the creation of a Langchain Chat object"""
import os
import glob
import langchain
from langchain.embeddings.openai import OpenAIEmbeddings

# from langchain.text_splitter import CharacterTextSplitter
from langchain.embeddings.mosaicml import MosaicMLInstructorEmbeddings
from langchain.chat_models import ChatOpenAI
from langchain.text_splitter import RecursiveCharacterTextSplitter
from langchain.chains import ConversationalRetrievalChain
from langchain.document_loaders import NotebookLoader
from langchain.vectorstores import FAISS
from langchain.document_loaders import TextLoader, DirectoryLoader
from openad.helpers.output import output_error, output_warning


## Creds clas for Watson X disabled currently
# class my_creds:
#    "Chat credentials Object"
#    DEFAULT_API = "https://workbench-api.res.ibm.com/v1"
#    api_key=None
#    api_endpoint=None
#    def __init__(
#        self,
#        api_key: str,
#        api_endpoint: str = DEFAULT_API,
#        ):
#        """
#        Instansiate the credentials object
#
#        Args:
#            api_key (str): The GENAI API Key
#            api_endpoint (str, optional): GENAI API Endpoint. Defaults to DEFAULT_API.
#        """
#        if api_key is None:
#            raise ValueError("api_key must be provided")
#        self.api_key = api_key
#        if api_endpoint is None:
#            raise ValueError("api_endpoint must be provided")
#        self.api_endpoint = api_endpoint


class Chatobject:
    """This is the Chat Object that is instantiated once per session"""

    llm_service = None
    llm_model = None
    organisation = None
    target = None
    API_key = None
    vector_db = None
    db_handle = None
    chat_history = []
    db_dir = "~/.vector_embed"
    document_folders = ["./"]
    document_types = ["**/*.txt", "**/*.ipynb", "**/*.run", "**/*.cdoc"]
    chat_template = """
    You are the Tell Me assistant that responds in a Helpful manner with responses like a Technical Documentation Writer.
    Here is the users current Request.
    ####

    {Question}

    ####

    Please consider the following Chat History as possible Context, If the same question is repeated see it as a request for more verbosity.
    ####
    {chat_history}
    ####
    """

    def __init__(
        self,
        target="OPENAPI",
        organisation=None,
        API_key=None,
        vector_db="FAISS",
        document_folders=["./"],
        document_types=document_types,
        db_dir_override=None,
        refresh_vector=False,
        llm_model="gpt-4",
        llm_service="OPENAI",
    ):
        self.organisation = organisation
        self.target = target
        self.API_key = API_key
        self.vector_db = vector_db
        self.db_handle = None
        self.document_folders = document_folders
        self.document_types = document_types
        self.llm_service = llm_service
        self.llm_model = llm_model

        langchain.embeddings.openai.api_key = self.organisation
        langchain.embeddings.organization = self.target

        if db_dir_override is not None and os.path.exists(os.path.expanduser(db_dir_override)) is True:
            self.db_dir = db_dir_override
        if self.vector_db == "FAISS":
            try:
                self.db_handle = self.load_faiss_db(refresh_vector)
            except Exception as e:  # pylint: disable=broad-exception-caught
                raise Exception(
                    f"the vector db {self.db_handle} Was not able to be loaded"
                ) from e  # pylint: disable=broad-exception-raised
        else:
            raise Exception(
                f"the vector db {self.db_handle} is not currently supported"
            )  # pylint: disable=broad-exception-raised

    def prime_chat_history(self, primer: str):
        """ "add the prompt tuning text primer to the chat"""
        self.chat_history.append((primer, "ok"))

    def load_faiss_db(self, refresh=True):
        """Load the Faiss Database Embeddings"""
        ###########################################################################
        # validation Testing
        main_db = None

        # Trap and return False if unable to generate emdeddings
        if self.llm_service == "OPENAI":
            try:
                embeddings = OpenAIEmbeddings(openai_api_key=self.API_key)
            except Exception as e:
                raise Exception(
                    "Error: cannot initialise embeddings, check API Key"
                ) from e  # pylint: disable=broad-exception-raised
        elif self.llm_service == "WATSONX":
            try:
                embeddings = MosaicMLInstructorEmbeddings(
                    endpoint_url=self.organisation, mosaicml_api_token=self.API_key
                )
            except Exception as e:
                raise Exception(
                    "Error: cannot initialise embeddings, check API Key"
                ) from e  # pylint: disable=broad-exception-raised
        # If not refreshing the database, check to see if the database exists
        if refresh is not True:
            try:
                if self.vector_db == "FAISS":
                    main_db = FAISS.load_local(
                        os.path.expanduser(self.db_dir + "/faiss_index"), embeddings
                    )  # pylint: disable=no-member
                return main_db
            except:  # pylint: disable=bare-except
                # if datatabase not there force a refresh
                refresh = True

        docs = []
        # Instruct the user as the tool has detected a change in underlying toolkt or workspace that it will update the FAISS index
        output_warning("Updating Embeddings for current Toolkits and Workspaces", return_val=False)
        try:
            # excluded_files=[]
            for i in self.document_folders:
                print("Embedding:", i)
                for j in self.document_types:
                    if j == "**/*.ipynb":
                        for file in glob.glob(i + "/*.ipynb"):
                            loader = NotebookLoader(
                                file,
                                include_outputs=False,
                                max_output_length=20,
                                remove_newline=False,
                            )
                            try:
                                documents = loader.load()
                                text_splitter = RecursiveCharacterTextSplitter(
                                    chunk_size=400, chunk_overlap=0, separators=[","]
                                )
                                docs.extend(text_splitter.split_documents(documents))
                            except:  # pylint: disable=bare-except
                                # excluded_files.append(file)
                                # Some notebook files are just not processable rather than notifying as user could have many,
                                # we skip over ones that cannot be processed
                                pass
                    elif j == "**/*.cdoc":
                        loader = DirectoryLoader(i, glob=j, loader_cls=TextLoader)
                        documents = loader.load()
                        text_splitter = RecursiveCharacterTextSplitter(
                            chunk_size=2000, chunk_overlap=100, separators=["\@"], keep_separator=False
                        )
                        docs.extend(text_splitter.split_documents(documents))
                    else:
                        loader = DirectoryLoader(i, glob=j, loader_cls=TextLoader)
                        documents = loader.load()
                        text_splitter = RecursiveCharacterTextSplitter(
                            chunk_size=700, chunk_overlap=100, separators=["\n"]
                        )
                        docs.extend(text_splitter.split_documents(documents))

            main_db = FAISS.from_documents(docs, embeddings)  # pylint: disable=no-member
            main_db.save_local(os.path.expanduser(self.db_dir + "/faiss_index"))

        except Exception as e:  # pylint: disable=broad-exception-caught
            output_error("Error in creating vector database " + str(e), return_val=False)
            return False
        return main_db

    def how_to_search(self, search: str):
        """Executing the Tell Me Function"""
        retriever = self.db_handle.as_retriever()
        if self.llm_service == "OPENAI":
            try:
                model = ChatOpenAI(
                    model_name=self.llm_model,
                    openai_api_key=self.API_key,
                )  # Other options 'ada' 'gpt-3.5-turbo' 'gpt-4'
            except Exception as e:  # pylint: disable=broad-exception-caught
                return output_error("Error Loading OPENAI Model see error Messsage : \n" + e, return_val=True)
        # this code exists to continue on once Watson has large enough model
        #   elif  self.llm_service == 'WATSONX':
        #       try: ### currently under Development
        #           sys.stdout = open(os.devnull, "w")
        #           sys.stderr = open(os.devnull, "w")
        #            # This client spins out "DEPRECATION" Warnings like nothing else so I am silencing it.
        #            from genai.credentials import Credentials
        #            from genai.extensions.langchain import LangChainInterface
        #            from genai.model import Model
        #            from genai.schemas import  GenerateParams
        #            creds = Credentials(api_key=self.API_key, api_endpoint=self.organisation)
        #            params = GenerateParams(decoding_method="greedy",max_new_tokens=None)
        #            model = LangChainInterface(model=self.llm_model,params=params  ,credentials=creds,verbose=True)
        #            sys.stdout = sys.__stdout__
        #            sys.stderr = sys.__stderr__
        #        except Exception as e:
        #            return  "Error Loading Watson X Model see error Messsage : \n"+e
        #      """

        # below code allows for multiple staged questions to be asked however this is not enabled in current version
        try:
            qa = ConversationalRetrievalChain.from_llm(model, retriever=retriever)
            questions = [search]
            answers = None
            for question in questions:
                try:
                    result = qa({"question": question, "chat_history": self.chat_history})
                except Exception as e:  # pylint: disable=broad-exception-caught
                    return output_error("Unable to Execute Request: " + str(e), return_val=True)
        except Exception as e:  # pylint: disable=broad-exception-caught
            return output_error("Failed Querying LLM " + str(e), return_val=True)
        try:
            self.chat_history.append((question, result["answer"]))
            if len(self.chat_history) > 3:
                try:
                    self.chat_history.remove(2)
                except Exception:  # pylint: disable=broad-exception-caught
                    pass
            answers = result["answer"]
        except Exception as e:  # pylint: disable=broad-exception-caught
            return output_error("Unable to Execute Request: " + str(e), return_val=True)
        return answers
