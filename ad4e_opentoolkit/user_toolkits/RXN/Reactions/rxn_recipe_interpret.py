

_tableformat = 'simple'
from ad4e_opentoolkit.helpers.output import output_text as output_text
from rxn4chemistry import RXN4ChemistryWrapper
list_of_reactions =[]

def get_include_lib(cmd_pointer):
    import importlib.util as ilu
    folder = cmd_pointer.toolkit_dir+'/RXN'+'/rxn_include.py'
    file = 'rxn_include'
    spec = ilu.spec_from_file_location(file, folder)
    rxn = ilu.module_from_spec(spec)
    spec.loader.exec_module(rxn)
    rxn_helper= rxn.rxn_helper
    return rxn_helper

from rdkit import Chem
from rdkit.Chem import AllChem
from typing import Dict, List


def interpret_recipe(inputs: dict, toolkit_dir, cmd_pointer):
   
    receipe = inputs["receipe"]
    
    rxn4chemistry_wrapper = cmd_pointer.login_settings['client'][cmd_pointer.login_settings['toolkits'].index('RXN') ]
    # Prepare the data query
  
   
    if cmd_pointer.notebook_mode == True:
        from IPython.display import display,Markdown
    try:
        return_result=[]
        
        actios_from_procedure_results = rxn4chemistry_wrapper.paragraph_to_actions(receipe)
        if cmd_pointer.notebook_mode == True:
            return_result.append( "***See the following actions from the Receipe:***\n")
        else:
            return_result.append("See the following actions from the Receipe:\n")

        for index, action in enumerate(actios_from_procedure_results['actions'], 1):
            
            if cmd_pointer.notebook_mode == True:
                return_result.append(f'{index}. {action}\n')
            else:
                return_result.append(f'{index}. {action}\n')
        if cmd_pointer.notebook_mode==True:
            return Markdown(''.join(return_result))
        else:
            return output_text('\n'+''.join(return_result),cmd_pointer=cmd_pointer)
    except  BaseException as e:
            raise BaseException("unable to to turn paragraph to actions:" +str(e)    )
            
