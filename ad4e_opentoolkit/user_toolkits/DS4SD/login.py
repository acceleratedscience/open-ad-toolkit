import os
from pathlib import Path
import deepsearch as ds
from deepsearch.cps.client.components.elastic import ElasticDataCollectionSource
from deepsearch.cps.queries import DataQuery
from deepsearch.cps.client.components.queries import RunQueryError
from ad4e_opentoolkit.helpers.output import msg, output_text, output_error, output_warning

config_blank = {"host": "None", "auth": {"username": "None", "api_key": "None"}, "verify_ssl": "false"}

# login

# https://cps.foc-deepsearch.zurich.ibm.com/projects/1234567890abcdefghijklmnopqrstvwyz123456/library


# Initialize the Deep Search client from the config file
# Input parameters for the example flow
def reset(cmd_pointer):
    if  os.path.isfile(os.path.expanduser(cmd_pointer.home_dir) + "/ds-auth.ext-v2.json"):
        os.remove(os.path.expanduser(cmd_pointer.home_dir) + "/ds-auth.ext-v2.json")

def login(cmd_pointer):
    # {'toolkits':[],'toolkits_details':{},'toolkits_api':[]}
    
    import time
    if not os.path.isfile(os.path.expanduser(cmd_pointer.home_dir) + "/ds-auth.ext-v2.json"):
        login_reset = True
    else:
        login_reset= False

        
    if 'DS4SD' not in cmd_pointer.login_settings['toolkits']:
        cmd_pointer.login_settings['toolkits'].append('DS4SD')
        cmd_pointer.login_settings['toolkits_details'].append({"type": "config_file", "session": "handle"})
        cmd_pointer.login_settings['toolkits_api'].append(None)
        cmd_pointer.login_settings['client'].append(None)
        cmd_pointer.login_settings['expiry'].append(None)
        x = cmd_pointer.login_settings['toolkits'].index('DS4SD')
    elif login_reset==False:
        import datetime
        from datetime import datetime, timezone
        now = datetime.now(timezone.utc)
        x = cmd_pointer.login_settings['toolkits'].index('DS4SD')
        now = datetime.timestamp(now)
        expiry_time = cmd_pointer.login_settings['expiry'][x]

        if expiry_time is not None and expiry_time > now:
            expiry_datetime = time.strftime('%a %b %e, %G  at %R', time.localtime(expiry_time))
            return True, expiry_datetime

    if not os.path.isfile(os.path.expanduser(cmd_pointer.home_dir) + "/ds-auth.ext-v2.json"):
        output_warning('Setting Authentication Details for Deep Search:',cmd_pointer)
        import readline
        if cmd_pointer.notebook_mode == False:
            import readline          
            output_text("Enter the Hostname: if the hostname is left blank it will default to 'https://sds.app.accelerate.science/' ",cmd_pointer=cmd_pointer)
            config_blank['host'] = input("\u001b[33mHostname: \u001b[0m")
            if config_blank['host'].strip()=='':
                config_blank['host']='https://sds.app.accelerate.science/'
            
            readline.remove_history_item(readline.get_current_history_length() - 1)
            config_blank['auth']['username'] = input("\u001b[33mEmail: \u001b[0m")
            readline.remove_history_item(readline.get_current_history_length() - 1)
            config_blank['auth']['api_key'] = input("\u001b[33mApi_key: \u001b[0m")
            readline.remove_history_item(readline.get_current_history_length() - 1)
        else:    
            output_text("Enter the Hostname: if the hostname is left blank it will default to 'https://sds.app.accelerate.science/' ",cmd_pointer=cmd_pointer)
            config_blank['host'] = input("Hostname: ")
            if config_blank['host'].strip()=='':
                config_blank['host']='https://sds.app.accelerate.science/'
            readline.remove_history_item(readline.get_current_history_length() - 1)
            config_blank['auth']['username'] = input("Email: ")
            readline.remove_history_item(readline.get_current_history_length() - 1)
            config_blank['auth']['api_key'] = input("Api_key: ")
            readline.remove_history_item(readline.get_current_history_length() - 1)
        import json
        with open(os.path.expanduser(cmd_pointer.home_dir) + "/ds-auth.ext-v2.json", 'w') as handle:
            json.dump(config_blank, handle)
            handle.close()
        output_text('config file generated.',cmd_pointer=cmd_pointer)
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

        import jwt
        import time

        # decoded_token = jwt.decode(bearer, verify=False,algorithms="HS256")
        decoded_token = jwt.decode(bearer, options={'verify_at_hash': False, 'verify_signature': False}, verify=False)

        # Extract expiry time from token payload
        expiry_time = decoded_token['exp']

        # Convert expiry time to a human-readable format
        expiry_datetime = time.strftime('%a %b %e, %G  at %R', time.localtime(expiry_time))
        cmd_pointer.login_settings['expiry'][x] = expiry_time

        return True, expiry_datetime
    except BaseException as err:
        output_error(msg('err_login', 'DS4SD',"Unable to connect to "+config.host, split=True), cmd_pointer)
        return False, None
