""" Interface libraries for working with the LLM API"""
import os
import re
import shutil
import glob
from openad.llm_assist.prime_chat import Chatobject
from openad.helpers.output import output_text, output_error
from openad.helpers.output import output_warning, output_success
from openad.app.global_var_lib import _repo_dir
from openad.app.global_var_lib import _meta_dir
from openad.helpers.credentials import load_credentials
from openad.helpers.credentials import write_credentials, get_credentials

# Constants
TRAINING_LLM_DIR = "/prompt_train/"
SUPPORTED_LLMS = ["WATSONX", "OPENAI"]
PROMPT_DIR = "~/.chat_embedded"
STANDARD_FILE_TYPES_EMBED = ["*.txt", "*.ipynb", "*.run", "*.cdoc"]
EXTENDED_FILE_TYPES_EMBED = ["**/*.txt", "**/*.ipynb", "**/*.run", "**/*.cdoc"]
NOTEBOOKS_DIR = "/../notebooks"
DEFAULT_SOURCES_LIST = []  #
# CHAT_PRIMER_old = """  In formatting the answer use markdown formatting syntax to highlight \
#      instances of commands, Parameters, Examples, Command Options and Command Clauses. Here is a correct version of a example code \
#       ``` search collection 'pubchem' for 'Ibuprofen' show (data)  ```. \
#       Only format one line at a time. No "\n" characters in codeblocks. Do not mention formatting in response. Always format the response.
#         Tell me """
CHAT_PRIMER = """  In answering the following Question:
            -Explain what any requested commands do
            -Provide All syntax, Options, Parameters and Examples when answering a question
            -Provide All Command Syntax, Clauses or Option Syntax in codeblock Markdown
            -Only format one line at a time. codeblocks should not go over lines
            -No "\n" characters in codeblocks
             Here is an correct version of an example command in codeblock format ``` search collection 'PubChem' for 'Ibuprofen' show (data)  ```. \
            Only Use the above as a guide to how to respond not as content for a response. 
            Always format the answer.
            Tell me """
CHAT_PRIMER_SUFFIX = """. Please provide information on options in response. Please format all code syntax, clauses and options in codeblock formatting on single lines in the response."""
CHAT_HISTORY_PRIMER = """When answering questions in the following chats,
             Answer like a technical helpful writer.
             When answering always show a Summary Description of any mentioned Command explaining what it does, codeblocked command syntax 
             and codeblocked command option Syntax as well as any parameters defined. format the codeblocks in markdown.
             Prioritise Help Text over Notebook examples. """


def create_train_repo(
    included_sources=DEFAULT_SOURCES_LIST, location_for_documents=PROMPT_DIR, document_types=STANDARD_FILE_TYPES_EMBED
):
    """Creates a Training Repository to build the embeddings from for the assistant"""
    try:
        if not os.path.exists(os.path.expanduser(location_for_documents)):
            os.mkdir(os.path.expanduser(location_for_documents))
        for file in glob.glob(os.path.expanduser(location_for_documents) + "/*"):
            os.remove(file)
    # catch any error that may be caused by OS in creating document localtion
    except Exception as err:  # pylint: disable=broad-exception-caught
        output_error(
            f"unable to create training repository directory {location_for_documents}" + str(err), return_val=False
        )
        return False

    for source in included_sources:
        for doc_type in document_types:
            for file in glob.glob(os.path.expanduser(source) + "/" + doc_type):
                shutil.copy(file, os.path.expanduser(location_for_documents) + "/")

    return True


# tell me / how do i Functional call


def how_do_i(cmd_pointer, parser):
    """Calling Tell Me Function"""
    cmd_pointer.settings["env_vars"]["llm_service"] = cmd_pointer.llm_service
    included_dirs = []
    if not os.path.exists(os.path.expanduser(PROMPT_DIR)):
        os.mkdir(os.path.expanduser(PROMPT_DIR))

    included_dirs.append(os.path.expanduser(_repo_dir + NOTEBOOKS_DIR))
    included_dirs.append(os.path.expanduser(_meta_dir + TRAINING_LLM_DIR))

    for name in cmd_pointer.settings["workspaces"]:
        included_dirs.append(f"{cmd_pointer.home_dir}/{name}")
    if cmd_pointer.llm_handle is None or cmd_pointer.refresh_vector is True:
        create_train_repo(included_sources=included_dirs, document_types=STANDARD_FILE_TYPES_EMBED)
        try:
            cmd_pointer.llm_handle = Chatobject(
                API_key=get_api_key(llm_name=cmd_pointer.llm_service, cmd_pointer=cmd_pointer)["auth"]["api_key"],
                organisation=get_api_key(llm_name=cmd_pointer.llm_service, cmd_pointer=cmd_pointer)["host"],
                document_folders=[os.path.expanduser(PROMPT_DIR)],
                document_types=EXTENDED_FILE_TYPES_EMBED,
                refresh_vector=cmd_pointer.refresh_vector,
                llm_service=cmd_pointer.llm_service,
                llm_model=cmd_pointer.llm_model,
            )
            if cmd_pointer.llm_handle is False:
                return False
        except Exception as e:  # pylint: disable=broad-exception-caught
            output_error("Problem Connecting to LLM: " + str(e), return_val=False, cmd_pointer=cmd_pointer)
            return False  # if there any other error in calling LLM handle e.g. network , loss of connection etc.

        cmd_pointer.refresh_vector = False
        cmd_pointer.settings["env_vars"]["refresh_help_ai"] = False
        # This puts Hostroy content that is instructional to the prompt on how to behave
        try:
            cmd_pointer.llm_handle.prime_chat_history(CHAT_HISTORY_PRIMER)
        except:
            #
            output_text(
                "Unable to Execute request. check LLM credentials and or Connectivity",
                return_val=False,
                cmd_pointer=cmd_pointer,
                pad=1,
                edge=True,
            )
            return False
    if cmd_pointer.notebook_mode is True:
        from halo import HaloNotebook as Halo  # pylint: disable=import-outside-toplevel
    else:
        from halo import Halo  # pylint: disable=import-outside-toplevel

    class Spinner(Halo):
        """custom spinner"""

        def __init__(self):
            # Alternative spinners:
            # simpleDotsScrolling, interval=100
            super().__init__(spinner="dots", color="white")

    newspin = Spinner()
    newspin.start("Processing Request ")
    # Now we are asking the prompt a Question

    try:
        text = cmd_pointer.llm_handle.how_to_search(CHAT_PRIMER + " ".join(parser["Chat_String"]) + CHAT_PRIMER_SUFFIX)
    except Exception as e:
        newspin.fail("Running Request Failed")
        output_text(
            "Unable to Execute request. check LLM credentials and or Connectivity",
            return_val=False,
            cmd_pointer=cmd_pointer,
            pad=1,
            edge=True,
        )
        return False
    newspin.succeed("See Answer Below.")
    text = clean_up_llm_text(cmd_pointer, text)

    if cmd_pointer.notebook_mode is True:
        return output_text(text, return_val=True, cmd_pointer=cmd_pointer, pad=1, edge=True)
    else:
        return output_text("\n" + text + "\n\n", return_val=True, cmd_pointer=cmd_pointer)


# sets the support llm model to use
def clean_up_llm_text(cmd_pointer, old_text):
    """This function cleans up text based on common LLM formatting and translates to our standard formatting"""

    text = old_text

    # LLM sometimes places the code type used inside the markdown section this simply removes it
    # Needs tidyup
    text = re.sub(r"\`\`\`python\n", r"```\n", text)
    text = re.sub(r"\`\`\`markdown\n", r"```\n", text)
    text = re.sub(r"\`markdown\n", r"`\n", text)
    text = re.sub(r"\`python\n", r"`\n", text)
    text = re.sub(r"\`\`\`plaintext\n", r"```\n", text)
    text = re.sub(r"\`plaintext\n", r"`\n", text)
    text = re.sub(r"\`\`\`python", r"```", text)
    text = re.sub(r"\`\`\`markdown", r"```", text)
    text = re.sub(r"\`markdown", r"`", text)
    text = re.sub(r"\`python", r"`", text)
    text = re.sub(r"\`\`\`plaintext", r"```", text)
    text = re.sub(r"\`plaintext", r"`", text)

    if cmd_pointer.notebook_mode is not True:
        # Needs optimising
        text = re.sub(r"\`\`\`\n[\s]+%openad ", r"```\n ", text)
        text = re.sub(r"\`\`\`\n\%openad ", r"```\n", text)
        text = re.sub(r"\`\`\`[\s]+%openad ", r"``` ", text)
        text = re.sub(r"\`\[\n\r\s]+%openad ", r"` ", text)
        text = re.sub(r"\`\%openad ", r"`", text)
        text = re.sub(r"\`[\s]+\%openad (.*?)\n", r"`\1\n", text)
        text = re.sub(r"\`[\s]+\%openad (.*?)", r"`\1", text)
        text = re.sub(r"[\s]+\%openad", r" ", text)

    # Replace ``` or ` with <cmd> bracing
    text = re.sub(r"\`\`\`(\n*?)(\s*?)(\%*?)([a-z]\n*[\s\S]*?)(\n*?)(\s*?)\`\`\`", r" <cmd>\3\4</cmd> ", text)
    text = re.sub(r"\`\`\`([a-z]*[\s\S]*?)\`\`\`", r" <cmd>\1</cmd> ", text)
    text = re.sub(r"\`([a-z]*[\s\S]*?)\`", r" <cmd>\1</cmd> ", text)
    text = re.sub(r"\`(\n*?)(\s*?)(\%*?)([a-z]\n*[\s\S]*?)(\n*?)(\s*?)\`", r" <cmd>\3\4</cmd> ", text)

    # nuance of llm instructued to use markdown

    if cmd_pointer.notebook_mode is not True:
        # LLMS Sometimes send back some interesting Markdown this is to translate it
        # It assumes that the LLM is not putting any formatting inside a code block
        text = re.sub(r"### (.*?) ###", r"<h2> \1 </h2>\n", text)
        text = re.sub(r"### (.*?)\n", r"<h2> \1 </h2>\n", text)
        text = re.sub(r"## (.*?) ##", r"<h2> \1 </h2>\n", text)
        text = re.sub(r"## (.*?)\n", r"<h2> \1 </h2>\n", text)
        text = re.sub(r"# (.*?)\n", r"<h1> \1 </h1>\n", text)
        # change bold to green
        text = re.sub(r"\*\*\*(.*?)\*\*\*", r"<green> \1 </green>", text)
        text = re.sub(r"\*\*(.*?)\*\*", r"<green> \1 </green>", text)
    return text


def set_llm(cmd_pointer, parser):
    """Set the current llm Model API"""
    llm_name = str(parser["llm_name"][0])
    if llm_name.upper() in SUPPORTED_LLMS:
        cmd_pointer.llm_service = llm_name.upper()
        cmd_pointer.settings["env_vars"]["llm_service"] = llm_name.upper()
        return output_success(
            " The Following has been set as the current llm service " + llm_name.upper(), cmd_pointer, pad=1
        )

    return output_text("The following is an invalid service " + llm_name.upper(), cmd_pointer, pad=1)


# removes the llm api key file


def clear_llm_auth(cmd_pointer, parser):  # pylint: disable=unused-argument
    """clears out the authentication file for the LLM"""
    if os.path.exists(f"{cmd_pointer.home_dir}/{cmd_pointer.llm_service.lower()}_api.cred"):
        os.remove(f"{cmd_pointer.home_dir}/{cmd_pointer.llm_service.lower()}_api.cred")
    return output_text("cleared API Auth", cmd_pointer, pad=1)


# retrieves the api key from the home directory of the  app


def get_api_key(llm_name, cmd_pointer):
    """get the nominated API key for the LLM"""
    api_config = load_credentials(f"{cmd_pointer.home_dir}/{llm_name.lower()}_api.cred")
    if api_config is None:
        output_warning("No Stored LLM Credentials", cmd_pointer=cmd_pointer, return_val=False)
        api_config = {"host": "None", "auth": {"username": "None", "api_key": "None"}, "verify_ssl": "false"}
        api_config = get_credentials(
            cmd_pointer=cmd_pointer, credentials=api_config, creds_to_set=["host", "auth:api_key"]
        )
        write_credentials(api_config, os.path.expanduser(cmd_pointer.home_dir + "/" + llm_name.lower() + "_api.cred"))
    return api_config
