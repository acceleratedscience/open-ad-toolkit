import os
from pathlib import Path
from rxn4chemistry import RXN4ChemistryWrapper
import json
config_blank = {"host": "None", "auth": {"username": "None", "api_key": "None"}, "verify_ssl": "false"}
from ad4e_opentoolkit.helpers.output import msg, output_text, output_error, output_warning
_default_url="https://rxn.app.accelerate.science"

# Initialize the rxn client from the config file
# Input parameters for the example flow


def get_include_lib(cmd_pointer):
    import importlib.util as ilu
    folder = cmd_pointer.toolkit_dir+'/RXN'+'/rxn_include.py'
    file = 'rxn_include'
    spec = ilu.spec_from_file_location(file, folder)
    rxn = ilu.module_from_spec(spec)
    spec.loader.exec_module(rxn)
    rxn_helper= rxn.rxn_helper()
    return rxn_helper

def reset(cmd_pointer):
    if  os.path.isfile(os.path.expanduser(cmd_pointer.home_dir) + "/rxn-auth.ext-v2.json"):
        os.remove(os.path.expanduser(cmd_pointer.home_dir) + "/rxn-auth.ext-v2.json")

def login(cmd_pointer):
  
    
    rxn_helper= get_include_lib(cmd_pointer)
    if not os.path.isfile(os.path.expanduser(cmd_pointer.home_dir) + "/rxn-auth.ext-v2.json"):
        login_reset = True
    else:
        login_reset= False
   
    expiry_time='No Expiry'

    if 'RXN' not in cmd_pointer.login_settings['toolkits']:
      
        cmd_pointer.login_settings['toolkits'].append('RXN')
        cmd_pointer.login_settings['toolkits_details'].append({"type": "config_file", "session": "handle"})
        cmd_pointer.login_settings['toolkits_api'].append(None)
        cmd_pointer.login_settings['client'].append(None)
        cmd_pointer.login_settings['expiry'].append(None)
        x = cmd_pointer.login_settings['toolkits'].index('RXN')
        cmd_pointer.login_settings['session_vars'].append(rxn_helper._RXN_VARS_TEMPLATE)
    
    elif login_reset == False:
        import datetime
        from datetime import datetime, timezone
        now = datetime.now(timezone.utc)
        x = cmd_pointer.login_settings['toolkits'].index('RXN')
        client = cmd_pointer.login_settings['client'][cmd_pointer.login_settings['toolkits'].index('RXN') ]
        try:
            print("\n")
            email = client.current_user()['response']["payload"]["email"]
            output_text("loggining in as "+email,cmd_pointer=cmd_pointer)
        except:
            output_error(msg('err_login', 'RXN',"Unable to connect to RXN Server", split=True), cmd_pointer)
            return False, None
        name,id=rxn_helper.get_current_project(cmd_pointer)
       
        if name!=cmd_pointer.settings['workspace']:
            rxn_helper.sync_up_workspace_name(cmd_pointer)
                            
        
        now = datetime.timestamp(now)
        return True, expiry_time
        

    if not os.path.isfile(os.path.expanduser(cmd_pointer.home_dir) + "/rxn-auth.ext-v2.json"):
        output_warning('Setting Authentication Details for RXN:',cmd_pointer)
        import readline
        print("\n")
        output_text("Enter the Hostname: if the hostname is left blank it will default to '" + _default_url + "' ", cmd_pointer=cmd_pointer)
            
        if cmd_pointer.notebook_mode == False:
            
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
        print("\n")
        output_text('config file generated.',cmd_pointer=cmd_pointer)
        
    try:
        with open(os.path.expanduser(cmd_pointer.home_dir) + "/rxn-auth.ext-v2.json", 'r') as handle:
            CONFIG_FILE=json.load(handle)
            
            handle.close()
        x = cmd_pointer.login_settings['toolkits'].index('RXN')
       
        client=RXN4ChemistryWrapper(api_key= CONFIG_FILE['auth']['api_key'], base_url= CONFIG_FILE['host'])
        email = client.current_user()['response']['payload']['email']
        print("\n")
        output_text("loggining in as "+email,cmd_pointer=cmd_pointer)
        cmd_pointer.login_settings['toolkits_api'][x] = CONFIG_FILE['auth']['api_key'] 
        cmd_pointer.login_settings['client'][x] = client
        #name,id=rxn_helper.get_current_project(cmd_pointer)
        rxn_helper.sync_up_workspace_name(cmd_pointer,reset=True)
        return True, expiry_time
    except BaseException as err:
         
         output_error(msg('err_login', 'RXN',"Unable to connect to "+ CONFIG_FILE['host'] , split=True), cmd_pointer)
         return False, None

