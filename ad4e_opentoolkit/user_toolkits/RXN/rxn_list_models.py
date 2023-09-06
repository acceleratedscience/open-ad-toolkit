

_tableformat = 'simple'
from  ad4e_opentoolkit.helpers.output import output_table as output_table


def list_models(inputs: dict, toolkit_dir, cmd_pointer):

    api_key =  cmd_pointer.login_settings['toolkits_api'][cmd_pointer.login_settings['toolkits'].index('RXN') ]
    rxn4chemistry_wrapper = cmd_pointer.login_settings['client'][cmd_pointer.login_settings['toolkits'].index('RXN') ]
    # Prepare the data query
   
    source_list = []
    try:
        x = rxn4chemistry_wrapper.list_models()
    except Exception as e:
        raise BaseException("unable to load models :"+str(e))
        return "unable to load models :"+e
    results=[]
    results2=[]
    for i in x:
        comma=''
        outstr=[]
        for ii in x[i]:
            outstr.append(ii['name'])
            
        results.append(i)
        results2.append(outstr)
    res_dict={"Models":results,"versions":results2}
    import pandas as pd
    df = pd.DataFrame.from_dict(res_dict)
    df.style.hide(axis='index')
    if cmd_pointer.notebook_mode == True:
        from IPython.display import display, HTML
        return HTML(df.to_html(index=False))
   
    output_table(df, tablefmt=_tableformat, headers=['Models','Versions'])
    #from tabulate import tabulate
    #return '\n'+tabulate(df, tablefmt=_tableformat, headers=['Models','Versions'], showindex=False)+'\n'
    