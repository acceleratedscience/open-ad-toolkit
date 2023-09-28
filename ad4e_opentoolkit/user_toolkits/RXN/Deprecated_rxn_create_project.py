

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
def create_project(inputs: dict, toolkit_dir, cmd_pointer):
    rxn_helper = get_include_lib(cmd_pointer)()
    
   
    if rxn_helper.validate_project(cmd_pointer,inputs['project_name'])==False:
        rxn4chemistry_wrapper = cmd_pointer.login_settings['client'][cmd_pointer.login_settings['toolkits'].index('RXN') ]
        # Prepare the data query
        result=False
        try:
            rxn4chemistry_wrapper.create_project(inputs['project_name'])
            result=True
        except BaseException as e:
            result=False
            raise BaseException("Unable to create project :"+str(e))
        text="Unable to create Project, try ```rxn list projects``` to select a available projects "
        if result ==False:
            if cmd_pointer.notebook_mode == True:
                    from IPython.display import display,Markdown
                    return rxn_helper.output("<fail>"+text+"</fail>",cmd_pointer)
            else:
                    print(text)
            return False
        else:
            try:
                text="Your RXN Project has successfully Created "+inputs['project_name']
                if cmd_pointer.notebook_mode == True:
                    from IPython.display import display,Markdown
                    return rxn_helper.output("<success>"+text+"</success>",cmd_pointer)
                else:
                    return text
            except BaseException as e:
                print(e)
                text="Unable to create Project, try 'rxn list projects' to select a valid project "
                if cmd_pointer.notebook_mode == True:
                    from IPython.display import display,Markdown
                    return rxn_helper.output("<fail>"+text+"</fail>",cmd_pointer)
                else:
                    print(text)
                return False
        return True
    else:
        
        text="Unable to create Project, name '"+inputs['project_name'] +"' already exists, try 'rxn list projects' to select a unique project name "
        if cmd_pointer.notebook_mode == True:
            from IPython.display import display,Markdown
            return rxn_helper.output("<fail>"+text+"</fail>",cmd_pointer)
        else:
            print(text)
        return False
