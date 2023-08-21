

_tableformat = 'simple'

from rxn4chemistry import RXN4ChemistryWrapper


# login

# https://cps.foc-deepsearch.zurich.ibm.com/projects/1234567890abcdefghijklmnopqrstvwyz123456/library

def list_projects(inputs: dict, toolkit_dir, cmd_pointer):

    api_key =  cmd_pointer.login_settings['toolkits_api'][cmd_pointer.login_settings['toolkits'].index('RXN') ]
    rxn4chemistry_wrapper = cmd_pointer.login_settings['client'][cmd_pointer.login_settings['toolkits'].index('RXN') ]
    # Prepare the data query
    source_list = []
    try:
        x = rxn4chemistry_wrapper.list_all_projects()['response']['payload']['content']
    except BaseException as e:
        raise BaseException(" Unable to contact RXN Server to list projects: "+str(e))
    import pandas as pd
    df = pd.DataFrame(x)
    #df['createdOn'] = pd.to_datetime(df['createdOn'], format='mixed')
    #df['modifiedOn'] = pd.to_datetime(df['modifiedOn'], format='mixed')
    df= df[['name','description','attempts']]
    
    
    df.style.hide(axis='index')
    if cmd_pointer.notebook_mode == True:
        from IPython.display import display, HTML
        return HTML(df.to_html(index=False))
        
    from tabulate import tabulate
    return '\n'+tabulate(df, tablefmt=_tableformat, headers="keys", showindex=False)+'\n'
    
    
   