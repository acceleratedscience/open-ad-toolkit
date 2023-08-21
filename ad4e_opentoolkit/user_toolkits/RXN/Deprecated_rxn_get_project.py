

_tableformat = 'simple'

from rxn4chemistry import RXN4ChemistryWrapper

from importlib.machinery import SourceFileLoader    

def get_include_lib(cmd_pointer):
    import importlib.util as ilu
    folder = cmd_pointer.toolkit_dir+'/RXN'+'/rxn_include.py'
    file = 'rxn_include'
    spec = ilu.spec_from_file_location(file, folder)
    rxn = ilu.module_from_spec(spec)
    spec.loader.exec_module(rxn)
    rxn_helper= rxn.rxn_helper
    return rxn_helper

# the sets the current RXN Project 
def get_project(inputs: dict, toolkit_dir, cmd_pointer):
    rxn_helper = get_include_lib(cmd_pointer)()

    rxn4chemistry_wrapper = cmd_pointer.login_settings['client'][cmd_pointer.login_settings['toolkits'].index('RXN') ]
    # Prepare the data query
   
    result,id,description =rxn_helper.get_current_project(cmd_pointer)
    
         
   
    if result ==None:
        text="No Projects have been set try ```rxn list projects``` to select a valid project "
        if cmd_pointer.notebook_mode == True:
                from IPython.display import display,Markdown
                return rxn_helper.output("<fail>"+text+"</fail>",cmd_pointer)
        else:
                print(text)
                return False
    else:
         
            if cmd_pointer.notebook_mode == True:
                text="Your RXN Project is set to <success>"+result+"</success>"
                from IPython.display import display,Markdown
                return rxn_helper.output(text,cmd_pointer)
            else:
                text="Your RXN Project has successfully been set to "+result
                return text
            
 


