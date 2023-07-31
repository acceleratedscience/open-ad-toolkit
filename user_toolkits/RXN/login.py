import os
from pathlib import Path
from rxn4chemistry import RXN4ChemistryWrapper
import json
config_blank = {"host": "None", "auth": {"username": "None", "api_key": "None"}, "verify_ssl": "false"}

# login

# https://cps.foc-deepsearch.zurich.ibm.com/projects/1234567890abcdefghijklmnopqrstvwyz123456/library


# Initialize the Deep Search client from the config file
# Input parameters for the example flow

import imp
import jwt
import time

def get_include_lib(cmd_pointer):
    import importlib.util as ilu
    folder = cmd_pointer.toolkit_dir+'/RXN'+'/rxn_include.py'
    file = 'rxn_include'
    spec = ilu.spec_from_file_location(file, folder)
    rxn = ilu.module_from_spec(spec)
    spec.loader.exec_module(rxn)
    rxn_helper= rxn.rxn_helper
    return rxn_helper

    


def login(cmd_pointer):
  
    #exec(open(cmd_pointer.toolkit_dir+'/RXN/rxn_include.py').read())
    #rxn_helper = importlib.import_module(cmd_pointer.toolkit_dir+'/RXN/rxn_include', "rxn_helper")
    
    rxn_helper= get_include_lib(cmd_pointer)
    
   
    expiry_time=None

    if 'RXN' not in cmd_pointer.login_settings['toolkits']:
      
        cmd_pointer.login_settings['toolkits'].append('RXN')
        cmd_pointer.login_settings['toolkits_details'].append({"type": "config_file", "session": "handle"})
        cmd_pointer.login_settings['toolkits_api'].append(None)
        cmd_pointer.login_settings['client'].append(None)
        cmd_pointer.login_settings['expiry'].append(None)
        x = cmd_pointer.login_settings['toolkits'].index('RXN')
        
        cmd_pointer.login_settings['session_vars'].append(rxn_helper._RXN_VARS_TEMPLATE)
    
    else:
        import datetime
        from datetime import datetime, timezone
        now = datetime.now(timezone.utc)
        x = cmd_pointer.login_settings['toolkits'].index('RXN')
     
        now = datetime.timestamp(now)
        return True, expiry_time
        #expiry_time = cmd_pointer.login_settings['expiry'][x - 1]

        #if expiry_time is not None and expiry_time > now:
        #    expiry_datetime = time.strftime('%a %b %e, %G  at %R', time.localtime(expiry_time))
        #    return True, expiry_datetime

    if not os.path.isfile(os.path.expanduser("~/.adccl") + "/rxn-auth.ext-v2.json"):

        if cmd_pointer.notebook_mode == False:
            print('\n'.join((
                "\n\u001b[31mConfiguration file for RXN not found:\u001b[0m",
                os.path.expanduser("~/.adccl") + "/rxn-auth.ext-v2.json",
                "\nPlease enter the following details to generate a new config file.\n"
            )))
            import readline
            config_blank['host'] = input("\u001b[33mHostname: \u001b[0m")
            readline.remove_history_item(readline.get_current_history_length() - 1)
            config_blank['auth']['api_key'] = input("\u001b[33mApi_key: \u001b[0m")
            readline.remove_history_item(readline.get_current_history_length() - 1)
        else:
            from IPython.display import display,Markdown
            display(Markdown(('<br>'.join((
                "<br> ***Configuration file for RXN not found:***",
                os.path.expanduser("~/.adccl") + "/rxn-auth.ext-v2.json",
                "<br> ***Please enter the following details to generate a new config file.***"
            )))))
            import readline
            config_blank['host'] = input("Hostname: ")
            readline.remove_history_item(readline.get_current_history_length() - 1)
            config_blank['auth']['api_key'] = input("Api_key: ")
            readline.remove_history_item(readline.get_current_history_length() - 1)


        
        with open(os.path.expanduser("~/.adccl") + "/rxn-auth.ext-v2.json", 'w') as handle:
            json.dump(config_blank, handle)
            handle.close()
        print('config file generated.')
    try:
        with open(os.path.expanduser("~/.adccl") + "/rxn-auth.ext-v2.json", 'r') as handle:
            CONFIG_FILE=json.load(handle)
            handle.close()
        
     
        client=RXN4ChemistryWrapper(api_key= CONFIG_FILE['auth']['api_key'], base_url= CONFIG_FILE['host'])
        
        client.list_all_projects()
        cmd_pointer.login_settings['toolkits_api'][x] = config_blank['auth']['api_key']
        cmd_pointer.login_settings['client'][x] = client
        
        return True, expiry_time
    except BaseException as err:
        print(err)
        pass
        # output_error(msg('err_login', 'DS4SD', err, split=True), cmd_pointer) # Don't have access to output_error here
        # print('\n'.join((
        #     '\n'
        #     '\u001b[31mWe were unable to log you into Deep Search.\u001b[0m',
        #     'Please check your credentials in your config file:',
        #     os.path.expanduser("~/.adccl") + "/ds-auth.ext-v2.json\n",
        #     f'\u001b[90mException: {err}\u001b[0m',
        #     ''
        # )))

        return False, None

