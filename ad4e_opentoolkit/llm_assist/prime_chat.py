#prime_chat

import os,sys,glob

import langchain 
from langchain.embeddings.openai import OpenAIEmbeddings
from langchain.text_splitter import CharacterTextSplitter
from langchain.text_splitter import RecursiveCharacterTextSplitter
from langchain.document_loaders import NotebookLoader
from langchain.vectorstores import FAISS
from langchain.document_loaders import TextLoader,DirectoryLoader
class my_creds:
        DEFAULT_API = "https://workbench-api.res.ibm.com/v1"
        api_key=None
        api_endpoint=None
        def __init__(
            self,
            api_key: str,
            api_endpoint: str = DEFAULT_API,
            ):
            """
            Instansiate the credentials object

            Args:
                api_key (str): The GENAI API Key
                api_endpoint (str, optional): GENAI API Endpoint. Defaults to DEFAULT_API.
            """
            if api_key is None:
                raise ValueError("api_key must be provided")
            self.api_key = api_key
            if api_endpoint is None:
                raise ValueError("api_endpoint must be provided")
            self.api_endpoint = api_endpoint

class  chat_object( ):
    import langchain
    llm_service     =   None
    llm_model       =   None
    organisation    =   None
    target          =   None
    API_key         =   None    
    vector_db       =   None
    db_handle       =   None
    chat_history    =   []
    db_dir          =   "~/.vector_embed"
    document_folders=   ["./"]
    document_types  =["**/*.txt","**/*.ipynb","**/*.run","**/*.cdoc"]

    
    def __init__(self,target='OPENAPI',organisation='org-V3VSRAXasFUnufPII8o1DIIk',API_key=None,vector_db='FAISS',document_folders=["./"],document_types=document_types,db_dir_override=None,refresh_vector=False,llm_model='gpt-4',llm_service='OPENAI'):
        self.organisation       =   organisation
        self.target             =   target
        self.API_key            =   API_key   
        self.vector_db          =   vector_db
        self.db_handle          =   None
        self.document_folders   = document_folders
        self.document_types     = document_types
        self.llm_service        =   llm_service 
        self.llm_model          =   llm_model  

        langchain.embeddings.openai.api_key=self.organisation
        langchain.embeddings.organization = self.target
        
        if db_dir_override!=None and os.path.exists(os.path.expanduser(db_dir_override)) == True:
             self.db_dir=db_dir_override    
       
        if self.vector_db == 'FAISS':
            try:
                self.db_handle  =  self.load_FAISS_db(refresh_vector)
                if self.db_handle==False:
                    
                    return False
            except  Exception as e:
                
                return False
        else:
            raise Exception(f"the vector db {self.db_handle} is not currently supported")
        
    def prime_chat_history(self,primer:str):
        self.chat_history.append((primer,"ok"))

    def load_FAISS_db(self,refresh=True):
        ###########################################################################
        #validation Testing
        delete_list=[]
        main_db=None

       
        # Trap and return False if unable to generate emdeddings
        if self.llm_service=='OPENAI' :
            try:
                embeddings = OpenAIEmbeddings(openai_api_key=self.API_key)
            except Exception as e:
                print(e)
                raise Exception("Error: cannot initialise embeddings, check API Key")
        elif self.llm_service=='WATSONX':
            try:
                from langchain.embeddings.mosaicml import MosaicMLInstructorEmbeddings
                embeddings = MosaicMLInstructorEmbeddings(endpoint_url=self.organisation ,mosaicml_api_token=self.API_key)
                
            except Exception as e:
                print(e)
                raise Exception("Error: cannot initialise embeddings, check API Key")
            
        
        # If not refreshing the database, check to see if the database exists 
        if refresh != True  : 
            try:
                if self.vector_db == 'FAISS':
                    main_db=FAISS.load_local(os.path.expanduser(self.db_dir+"/faiss_index"), embeddings)
                return main_db
            except:
                refresh = True
           
        for i in self.document_folders:
            if not os.path.exists(os.path.expanduser(i)):
                delete_list.append(i)
        for i in delete_list:
            self.document_folders.remove(i)
            print(f"Warning: removing folder '{i}' It does not exist or the application does not have access to it")

        
        docs=[]
        print("Loading Document Index:")
        
        try:
            for i in self.document_folders:
                for j in self.document_types:
                    if j =="**/*.ipynb":
                        for file in glob.glob(i+'/*.ipynb' ):
                            
                            
                            loader = NotebookLoader(file,include_outputs=False,max_output_length=20,remove_newline=False,)
                            
                            try:
                                documents = loader.load()
                            
                                #text_splitter = CharacterTextSplitter(chunk_size=3000, chunk_overlap=300,separators=[ "\n"])
                                text_splitter = RecursiveCharacterTextSplitter(chunk_size=400, chunk_overlap=0, separators=[","])
                                docs.extend(text_splitter.split_documents(documents))
                            except:
                                pass
                    elif j=="**/*.cdoc":
                        loader = DirectoryLoader(i,glob=j,loader_cls=TextLoader)
                        documents = loader.load()
                    #text_splitter = CharacterTextSplitter(chunk_size=3000, chunk_overlap=0)
                        #text_splitter = CharacterTextSplitter(chunk_size=3000, chunk_overlap=300,separators=[ "\n"])
                        text_splitter = RecursiveCharacterTextSplitter ( chunk_size=700, chunk_overlap=0, separators=["\@"])
                        
                        docs.extend(text_splitter.split_documents(documents))
                    
                    else:
                        loader = DirectoryLoader(i,glob=j,loader_cls=TextLoader)
                        documents = loader.load()
                    #text_splitter = CharacterTextSplitter(chunk_size=3000, chunk_overlap=0)
                        #text_splitter = CharacterTextSplitter(chunk_size=3000, chunk_overlap=300,separators=[ "\n"])
                        text_splitter = RecursiveCharacterTextSplitter(chunk_size=700, chunk_overlap=100, separators=["\n"])
                        docs.extend(text_splitter.split_documents(documents))
          
           
            main_db = FAISS.from_documents(docs, embeddings)
           
            main_db.save_local(os.path.expanduser(self.db_dir+"/faiss_index"))
           
        except Exception as e:
            print("error in creating vector database")
            print(e)
            return False
       
        
        return main_db
    

    def how_to_search(self,search:str):
        
 
        from langchain.chains import ConversationalRetrievalChain
        from langchain import PromptTemplate, LLMChain
        
        retriever=self.db_handle.as_retriever()
        
        if self.llm_service == 'OPENAI':
            try:

                from langchain.chat_models import ChatOpenAI
               
                model = ChatOpenAI(model_name=self.llm_model,openai_api_key=self.API_key)  # Other options 'ada' 'gpt-3.5-turbo' 'gpt-4',
                
            except Exception as e:
                
                return  "Error Loading OPENAI Model see error Messsage : \n"+e
        elif  self.llm_service == 'WATSONX':
            try:
                sys.stdout = open(os.devnull, "w")
                sys.stderr = open(os.devnull, "w")
                # This client spins out "DEPRECATION" Warnings like nothing else so I am silencing it.
                from genai.credentials import Credentials
                from genai.extensions.langchain import LangChainInterface
                from genai.model import Model
                from genai.schemas import  GenerateParams
                creds = Credentials(api_key=self.API_key, api_endpoint=self.organisation)
                params = GenerateParams(decoding_method="greedy",max_new_tokens=None)
                model = LangChainInterface(model=self.llm_model,params=params  ,credentials=creds)
                sys.stdout = sys.__stdout__
                sys.stderr = sys.__stderr__
            except Exception as e:
                return  "Error Loading Watson X Model see error Messsage : \n"+e
            
            
            
        try:    
            qa = ConversationalRetrievalChain.from_llm(model, retriever=retriever)
            questions=[search]
            answers=None
            

            for question in questions:
                
                try:
                    result = qa({"question": question, "chat_history": self.chat_history})
                except Exception as e:
                    print(e)
                    return 'Fail'
        except Exception as e:
            print(e)
            return 'Fail'

        try:            
            self.chat_history.append((question, result["answer"]))
            if len(self.chat_history) > 3:
                try:
                    self.chat_history.remove(2)
                except:
                    pass
            answers=result["answer"]
        except  Exception as e:
            print(e)
        return answers

