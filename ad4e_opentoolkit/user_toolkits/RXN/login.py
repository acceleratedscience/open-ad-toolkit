"""Login Procedure for the RXN Toolkit"""

import os
import json
import datetime
import importlib.util as ilu
from datetime import datetime, timezone
import readline
from ad4e_opentoolkit.helpers.output import msg, output_text, output_error, output_warning

from rxn4chemistry import RXN4ChemistryWrapper

config_blank = {"host": "None", "auth": {"username": "None", "api_key": "None"}, "verify_ssl": "false"}

_default_url="https://rxn.app.accelerate.science"

# Initialize the rxn client from the config file
# Input parameters for the example flow


def get_include_lib(cmd_pointer):
    """load rxn include library"""
    folder = cmd_pointer.toolkit_dir+'/RXN'+'/rxn_include.py'
    file = 'rxn_include'
    spec = ilu.spec_from_file_location(file, folder)
    rxn = ilu.module_from_spec(spec)
    spec.loader.exec_module(rxn)
    rxn_helper= rxn.rxn_helper()
    return rxn_helper

def reset(cmd_pointer):
    """remove on reset signal, the underlying API key file to trigger reset"""
    if  os.path.isfile(os.path.expanduser(cmd_pointer.home_dir) + "/rxn-auth.ext-v2.json"):
        os.remove(os.path.expanduser(cmd_pointer.home_dir) + "/rxn-auth.ext-v2.json")

def login(cmd_pointer):
    """logs onto the RXN service"""
    rxn_helper= get_include_lib(cmd_pointer)
    if not os.path.isfile(os.path.expanduser(cmd_pointer.home_dir) + "/rxn-auth.ext-v2.json"):
        login_reset = True
    else:
        login_reset= False
    first_login = False
   
    if 'RXN' not in cmd_pointer.login_settings['toolkits']: # If we have not already logged onto the RXN service in this Session
        first_login =True
        cmd_pointer.login_settings['toolkits'].append('RXN')
        cmd_pointer.login_settings['toolkits_details'].append({"type": "config_file", "session": "handle"})
        cmd_pointer.login_settings['toolkits_api'].append(None)
        cmd_pointer.login_settings['client'].append(None)
        cmd_pointer.login_settings['expiry'].append(None)
        x = cmd_pointer.login_settings['toolkits'].index('RXN')
        cmd_pointer.login_settings['session_vars'].append(rxn_helper._RXN_VARS_TEMPLATE) #pylint: disable=protected-access
    
    elif login_reset is False: # If a login reset has been issued or there is no authentication file
        
        now = datetime.now(timezone.utc)
        x = cmd_pointer.login_settings['toolkits'].index('RXN')
        client = cmd_pointer.login_settings['client'][cmd_pointer.login_settings['toolkits'].index('RXN') ]
        try:
            if login_reset is True or first_login is True:
                email = client.current_user()['response']["payload"]["email"]
                output_text("<success>loggining in as: </success> "+email,cmd_pointer=cmd_pointer, return_val=False)
        except Exception: #pylint: disable=broad-exception-caught
            output_error(msg('err_login', 'RXN',"Unable to connect to RXN Server", split=True),cmd_pointer=cmd_pointer, return_val=False)
            return False, None
        name,id=rxn_helper.get_current_project(cmd_pointer)
       
        if name!=cmd_pointer.settings['workspace']:
            rxn_helper.sync_up_workspace_name(cmd_pointer)
                            
        
        now = datetime.timestamp(now)
        return True,None
        
# if no Authentication file ask for authentication details and create
    if not os.path.isfile(os.path.expanduser(cmd_pointer.home_dir) + "/rxn-auth.ext-v2.json"): 
        output_warning('Setting Authentication Details for RXN:',cmd_pointer=cmd_pointer, return_val=False)
        output_text("Enter the Hostname: if the hostname is left blank it will default to '" + _default_url + "' ", cmd_pointer=cmd_pointer, return_val=False)
            
        if cmd_pointer.notebook_mode is False:
            config_blank['host'] = input("\u001b[33mHostname: \u001b[0m")
            if config_blank['host'].strip()=='':
                config_blank['host']=_default_url
            readline.remove_history_item(readline.get_current_history_length() - 1)
            config_blank['auth']['api_key'] = input("\u001b[33mApi_key: \u001b[0m")
            readline.remove_history_item(readline.get_current_history_length() - 1)
        else:
            config_blank['host']=input("Hostname: ")
            if config_blank['host'].strip()=='':
                config_blank['host']=_default_url   
            readline.remove_history_item(readline.get_current_history_length() - 1)
            config_blank['auth']['api_key'] = input("Api_key: ")
            readline.remove_history_item(readline.get_current_history_length() - 1)

        with open(os.path.expanduser(cmd_pointer.home_dir) + "/rxn-auth.ext-v2.json", 'w') as handle:
            json.dump(config_blank, handle)
            handle.close()
        output_text('config file generated.',cmd_pointer=cmd_pointer, return_val=False)
    try:
        with open(os.path.expanduser(cmd_pointer.home_dir) + "/rxn-auth.ext-v2.json", 'r') as handle:
            CONFIG_FILE=json.load(handle)
            
            handle.close()
        x = cmd_pointer.login_settings['toolkits'].index('RXN')
       
        client=RXN4ChemistryWrapper(api_key= CONFIG_FILE['auth']['api_key'], base_url= CONFIG_FILE['host'])
        email = client.current_user()['response']['payload']['email']
        if login_reset is True or first_login is True:
            output_text("<success>loggining in as:</success> "+email,cmd_pointer=cmd_pointer, return_val=False)
        cmd_pointer.login_settings['toolkits_api'][x] = CONFIG_FILE['auth']['api_key'] 
        cmd_pointer.login_settings['client'][x] = client 
        rxn_helper.sync_up_workspace_name(cmd_pointer,reset=True)
        return True, None
    except Exception: #pylint: disable=broad-exception-caught
        output_error(msg('err_login', 'RXN',"Unable to connect to "+ CONFIG_FILE['host'] , split=True),cmd_pointer=cmd_pointer, return_val=False)
        return False, None

