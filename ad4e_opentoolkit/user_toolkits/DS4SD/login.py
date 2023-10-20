import os
from  getpass import getpass
import datetime
from datetime import datetime, timezone
import time
from pathlib import Path
import requests
import jwt
import deepsearch as ds
from ad4e_opentoolkit.helpers.output import msg, output_text, output_error, output_warning

_default_url= 'https://sds.app.accelerate.science/'
config_blank = {"host": "None", "auth": {"username": "None", "api_key": "None"}, "verify_ssl": "false"}

# Initialize the Deep Search client from the config file
# Input parameters for the example flow
def reset(cmd_pointer):
    if  os.path.isfile(os.path.expanduser(cmd_pointer.home_dir) + "/ds-auth.ext-v2.json"):
        os.remove(os.path.expanduser(cmd_pointer.home_dir) + "/ds-auth.ext-v2.json")

def uri_valid(url:str)->bool:
    try:
        request = requests.get(url, stream=True,timeout=10)
    except:# pylint: disable=bare-except
        return False
    if request.status_code == 200:
        return True
    else:
        return False
def login(cmd_pointer):
    """Logs the Framework into the DS4SD service as defined by the user"""
    if not os.path.isfile(os.path.expanduser(cmd_pointer.home_dir) + "/ds-auth.ext-v2.json"):
        login_reset = True
    else:
        login_reset= False
    first_login =False
    if 'DS4SD' not in cmd_pointer.login_settings['toolkits']:
        cmd_pointer.login_settings['toolkits'].append('DS4SD')
        cmd_pointer.login_settings['toolkits_details'].append({"type": "config_file", "session": "handle"})
        cmd_pointer.login_settings['toolkits_api'].append(None)
        cmd_pointer.login_settings['client'].append(None)
        cmd_pointer.login_settings['expiry'].append(None)
        x = cmd_pointer.login_settings['toolkits'].index('DS4SD')
        first_login =True
    elif login_reset is False:
        now = datetime.now(timezone.utc)
        x = cmd_pointer.login_settings['toolkits'].index('DS4SD')
        now = datetime.timestamp(now)
        expiry_time = cmd_pointer.login_settings['expiry'][x]

        if expiry_time is not None and expiry_time > now:
            expiry_datetime = time.strftime('%a %b %e, %G  at %R', time.localtime(expiry_time))
            return True, expiry_datetime

    if not os.path.isfile(os.path.expanduser(cmd_pointer.home_dir) + "/ds-auth.ext-v2.json"):
        output_warning('Setting Authentication Details for Deep Search:', cmd_pointer)
        import readline
        if cmd_pointer.notebook_mode == False:
            import readline          
            output_text("Enter the URL / Hostname: if the hostname is left blank it will default to '"+_default_url+"' ", cmd_pointer=cmd_pointer,return_val=False)
            config_blank['host'] = input("\u001b[33m URL / Hostname: \u001b[0m")
            if config_blank['host'].strip()=='':
                config_blank['host']=_default_url
            elif uri_valid(config_blank['host']) ==False:
                output_error("Invalid URL Provided please check URL or VPN and try again",cmd_pointer=cmd_pointer,return_val=False)
                return False
            
            readline.remove_history_item(readline.get_current_history_length() - 1)
            config_blank['auth']['username'] = input("\u001b[33mEmail: \u001b[0m")
            readline.remove_history_item(readline.get_current_history_length() - 1)
            config_blank['auth']['api_key'] = getpass("\u001b[33mApi_key: \u001b[0m")
            readline.remove_history_item(readline.get_current_history_length() - 1)
        else:    
            output_text("Enter the URL / Hostname: if the hostname is left blank it will default to '" + _default_url + "' ", cmd_pointer=cmd_pointer,return_val=False) 
            config_blank['host'] = input("URL / Hostname: ")
            if config_blank['host'].strip()=='':
                config_blank['host']=_default_url
            if config_blank['host'].strip()=='':
                config_blank['host']=_default_url
            elif uri_valid(config_blank['host']) ==False:
                output_error("Invalid URL Provided please check URL or VPN and try again",cmd_pointer=cmd_pointer,return_val=False)
                return False
            readline.remove_history_item(readline.get_current_history_length() - 1)
            config_blank['auth']['username'] = input("Email: ")
            readline.remove_history_item(readline.get_current_history_length() - 1)
            config_blank['auth']['api_key'] = getpass("Api_key: ")
            readline.remove_history_item(readline.get_current_history_length() - 1)
        import json
        with open(os.path.expanduser(cmd_pointer.home_dir) + "/ds-auth.ext-v2.json", 'w') as handle:
            json.dump(config_blank, handle)
            handle.close()
        output_text('config file generated.',cmd_pointer=cmd_pointer,return_val=False)
    try:
        CONFIG_FILE = Path(os.path.expanduser(cmd_pointer.home_dir) + "/ds-auth.ext-v2.json")
        x = cmd_pointer.login_settings['toolkits'].index('DS4SD')
        config = ds.DeepSearchConfig.parse_file(CONFIG_FILE)

        client = ds.CpsApiClient(config)
        
        api = ds.CpsApi(client)

        cmd_pointer.login_settings['toolkits_api'][x ] = api
        cmd_pointer.login_settings['client'][x ] = client

        cb = client.bearer_token_auth
        bearer = cb.bearer_token

        # decoded_token = jwt.decode(bearer, verify=False,algorithms="HS256")
        decoded_token = jwt.decode(bearer, options={'verify_at_hash': False, 'verify_signature': False}, verify=False)

        # Extract expiry time from token payload
        expiry_time = decoded_token['exp']

        # Convert expiry time to a human-readable format
        expiry_datetime = time.strftime('%a %b %e, %G  at %R', time.localtime(expiry_time))
        if login_reset is True or first_login is True:
                output_text("<success>logged into DS4SD </success>",cmd_pointer=cmd_pointer, return_val=False)
        cmd_pointer.login_settings['expiry'][x] = expiry_time
        
        return True, expiry_datetime
    except Exception: # pylint: disable=broad-exception-caught
        output_error(msg('err_login', 'DS4SD',"Unable to connect to " + config.host, split=True), cmd_pointer=cmd_pointer,return_val=False)
        return False, None
