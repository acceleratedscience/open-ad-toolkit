from ad4e_opentoolkit.llm_assist.prime_chat import chat_object
import os
import shutil
import glob
from ad4e_opentoolkit.core.help import OpenadHelp
from ad4e_opentoolkit.helpers.output import msg, output_text, output_error, output_warning, output_success, output_table
from ad4e_opentoolkit.app.global_var_lib import _repo_dir
import pickle
import re
import readline
from ad4e_opentoolkit.app.global_var_lib import _meta_dir

_training_dir = '/prompt_train/'
_supported_llms = ['WATSONX', 'OPENAI']
_prompt_dir = '~/.chat_embedded'
_standard_types = ['*.txt', '*.ipynb', '*.run', '*.cdoc']
_extended_types = ['**/*.txt', '**/*.ipynb', '**/*.run', '**/*.cdoc']


def create_train_repo(included_sources=[], location_for_documents=_prompt_dir, document_types=_standard_types):
    try:

        if not os.path.exists(os.path.expanduser(location_for_documents)):
            os.mkdir(os.path.expanduser(location_for_documents))
    except BaseException:
        raise Exception(f"unable to create training repository directory {location_for_documents}")

    for file in glob.glob(os.path.expanduser(location_for_documents) + "/*"):

        os.remove(file)

    for source in included_sources:
        for type in document_types:
            for file in glob.glob(os.path.expanduser(source) + "/" + type):
                shutil.copy(file, os.path.expanduser(location_for_documents) + '/')

    return True
# tell me / how do i Functional call


def how_do_i(cmd_pointer, parser):
    cmd_pointer.settings['env_vars']['llm_service'] = cmd_pointer.llm_service

    included_dirs = [os.path.expanduser(_prompt_dir)]
    try:
        os.mkdir(os.path.expanduser(_prompt_dir))
    except BaseException:
        pass
    included_dirs.append(os.path.expanduser(_repo_dir + '/../notebooks'))
    included_dirs.append(os.path.expanduser(_meta_dir + _training_dir))
    for name in cmd_pointer.settings['workspaces']:
        included_dirs.append(f'{cmd_pointer.home_dir}/{name}')
    if cmd_pointer.llm_handle is None or cmd_pointer.refresh_vector == True:
        create_train_repo(included_sources=included_dirs, document_types=_standard_types)
        try:

            cmd_pointer.llm_handle = chat_object(API_key=get_api_key(cmd_pointer.llm_service, cmd_pointer)['auth']['api_key'],
                                                 organisation=get_api_key(cmd_pointer.llm_service, cmd_pointer)['host'], document_folders=[os.path.expanduser(_prompt_dir)], document_types=_extended_types, refresh_vector=cmd_pointer.refresh_vector, llm_service=cmd_pointer.llm_service, llm_model=cmd_pointer.llm_model)
            if cmd_pointer.llm_handle == False:
                return False
        except BaseException:

            return False

        cmd_pointer.refresh_vector = False
        cmd_pointer.settings['env_vars']['refresh_help_ai'] = False
        cmd_pointer.llm_handle.prime_chat_history('When answering questions in the following chats, Answer like a technical help writer \
        ,and always show the Command description and syntax including options. Note: %openad is the notebook magic command prompt. If an answer includes "%openad" note that these commands are for use in Notebooks ')
        # cmd_pointer.llm_handle.prime_chat_history('When answering questions in the following chats, Answer like a technical help writer Showing Syntax and examples. When answering always interpret  all Pyparsing_Command_Definitions using python pyparsing  and display only the matching user syntax without mentioning pyparsing at all, never mention pyparsing in answers\
        # , and always show the full command  then underneath bullet point syntax clauses highlighting required and optional syntax.')

    if cmd_pointer.notebook_mode == True:
        chat_primer = "Responding using Markdown format as if you are a Helpful Technical Writer, Tell me "

    else:
        chat_primer = 'Responding assuming a command Line output format as if you are a Helpful Technical Writer,  Tell Me '

    if cmd_pointer.notebook_mode == True:
        import IPython.display
        return IPython.display.Markdown(cmd_pointer.llm_handle.how_to_search(chat_primer + " ".join(parser['Chat_String'])))
    # interperate=' and in answer like a technical help writer and interperate all python pyparsing statements as Domain specific language should be entered by the user, never show the pyparsing syntax, and always show the full command and bullet point syntax clauses highlighting required and optional syntax'
    # return output_text( cmd_pointer.llm_handle.how_to_search("tell me " +" ".join(parser['Chat_String'])+interperate), cmd_pointer, pad=1, edge=True)

    # interperate=' and in answer like a technical help writer and interperate all python pyparsing statements as Domain specific language should be entered by the user, never show the pyparsing syntax, and always show the full command and bullet point syntax clauses highlighting required and optional syntax'
    # print(cmd_pointer.llm_handle.how_to_search(chat_primer +" ".join(parser['Chat_String'])))
    if not cmd_pointer.notebook_mode and not cmd_pointer.api_mode:
        text = re.sub(r'```(.*?)```', r'<cmd>\1</cmd>', cmd_pointer.llm_handle.how_to_search(chat_primer + " ".join(parser['Chat_String'])))
    else:
        text = cmd_pointer.llm_handle.how_to_search(chat_primer + " ".join(parser['Chat_String']))
    return output_text(text, return_val=True, cmd_pointer=cmd_pointer, pad=1, edge=True)

# sets the support llm model to use


def set_llm(cmd_pointer, parser):
    try:
        llm_name = str(parser['llm_name'][0])

        if llm_name.upper() in _supported_llms:
            cmd_pointer.llm_service = llm_name.upper()
            cmd_pointer.settings['env_vars']['llm_service'] = llm_name.upper()
            return output_text(' The Following has been set as the current llm service ' + llm_name.upper(), cmd_pointer, pad=1)
        else:
            return output_text('The following is an invalid service ' + llm_name.upper(), cmd_pointer, pad=1)
    except Exception as e:
        print(e)
        return output_text("error setting llm " + llm_name.upper(), cmd_pointer, pad=1)

# removes the llm api key file


def clear_llm_auth(cmd_pointer, parser):
    try:
        os.remove(os.path.expanduser(cmd_pointer.home_dir + "/" + cmd_pointer.llm_service.lower() + "_api.json"))

    except Exception as e:
        print(e)
        pass
    return output_text('cleared API Auth', cmd_pointer, pad=1)

# retrieves the api key from the home directory of the  app


def get_api_key(llm_name, cmd_pointer):
    try:
        api_config = load_api_registry(os.path.expanduser(cmd_pointer.home_dir + '/' + llm_name.lower() + "_api.json"))
    except Exception as e:
        print(e)
    if api_config is None:
        api_config = {"host": "None", "auth": {"username": "None", "api_key": "None"}, "verify_ssl": "false"}
        print('\n'.join((
            "\n\u001b[31m API Key Not found:\u001b[0m",
            os.path.expanduser(cmd_pointer.home_dir) + "/" + llm_name.lower() + "_api.json",
            "\nPlease enter the Organization: \n"
        )))

        api_config['host'] = input("\u001b[33m Hostname: \u001b[0m")
        readline.remove_history_item(readline.get_current_history_length() - 1)
        api_config['auth']['api_key'] = input("\u001b[33m Api_key: \u001b[0m")
        readline.remove_history_item(readline.get_current_history_length() - 1)
        write_api_registry(api_config, os.path.expanduser(cmd_pointer.home_dir + '/' + llm_name.lower() + "_api.json"))

    return api_config


# Load the user's api login data.
def load_api_registry(location):
    try:
        with open(location, 'rb') as handle:
            registry = pickle.loads(handle.read())
    except BaseException:
        registry = None
    return registry


# Dumps the LLM api details to home directory of app
def write_api_registry(registry: dict, location):

    try:
        with open(location, 'wb') as handle:
            settings = pickle.dump(registry, handle)
    except BaseException:
        return False
    return True
