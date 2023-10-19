

_tableformat = 'simple'

try:
    from rxn4chemistry import RXN4ChemistryWrapper
except:
    print("error loading rxn4chemistry")
    raise BaseException("error loading rxn4chemistry")
from time import sleep
def get_include_lib(cmd_pointer):
    import importlib.util as ilu
    folder = cmd_pointer.toolkit_dir+'/RXN'+'/rxn_include.py'
    file = 'rxn_include'
    spec = ilu.spec_from_file_location(file, folder)
    rxn = ilu.module_from_spec(spec)
    spec.loader.exec_module(rxn)
    rxn_helper= rxn.rxn_helper()
    return rxn_helper

# login
from rdkit import Chem
from rdkit.Chem import AllChem

def get_reaction_from_smiles(reaction_smiles: str) -> Chem.rdChemReactions.ChemicalReaction:
    return AllChem.ReactionFromSmarts(reaction_smiles, useSmiles=True)
# https://cps.foc-deepsearch.zurich.ibm.com/projects/1234567890abcdefghijklmnopqrstvwyz123456/library

def predict_reaction(inputs: dict, toolkit_dir, cmd_pointer):
    rxn_helper=get_include_lib(cmd_pointer)
    rxn_helper.sync_up_workspace_name(cmd_pointer)
    name,id = rxn_helper.get_current_project(cmd_pointer)
    if cmd_pointer.notebook_mode == True:
        import IPython
        from halo import HaloNotebook as Halo
    else:
        from halo import Halo

    if name == None and cmd_pointer.api_mode==False:
        if cmd_pointer.notebook_mode==True:
            from IPython.display import display, Markdown
            display(Markdown(" No current RXN project selected ,`set rxn project <project name>` to set your project before proceeding."))
            display(Markdown("Select from the Below Projects or create a new."))
            display(rxn_helper.get_all_projects(cmd_pointer)[['name','description']])
        else:
            print(" No current RXN project selected ,`set rxn project <project name>` to set your project before proceeding. ")
            print(" Select from the Below Projects or create a new.")
            print(rxn_helper.get_all_projects(cmd_pointer)[['name','description']])
        return False


    rxn4chemistry_wrapper = cmd_pointer.login_settings['client'][cmd_pointer.login_settings['toolkits'].index('RXN') ]
    # Prepare the data query
    rxn_helper = get_include_lib(cmd_pointer) 
    molecule = inputs["molecule"]
    val='val'
    predict_id=None
    ai_model='2020-08-10'
    if 'prediction_id' in inputs:
       predict_id = inputs["prediction_id"][val]
    if 'ai_model' in inputs:
        ai_model = inputs["ai_model"][val]
    if cmd_pointer.notebook_mode == True:
        from IPython.display import display,Markdown
    error_list=[]
    if 'use_saved'in inputs:
         use_saved = True
    else:
        use_saved=False
    
    
    for i in molecule.split('.'):
        if not rxn_helper.valid_smiles(str(i)):
            error_list.append(i)
        
    if len(error_list)>0:
        if cmd_pointer.notebook_mode == True:
            import pandas as pd
            df = pd.DataFrame(error_list,columns=['smiles'])
            
            display(Markdown("***Error:*** The following invalid were Smiles Supplied:"))
            display(df)
        else:
            print("Error: The following invalid were Smiles Supplied:")
            import pandas as pd
            df = pd.DataFrame(error_list,columns=['smiles'])
            print(df)
    ##################################################################################################################
    # check for cached reaction
    entry_2=[]

    for i in molecule.split('.'):
        entry_2.append(Chem.MolToSmiles(Chem.MolFromSmiles(i),canonical=True))
        dot='.'

    result = rxn_helper.retrieve_cache(cmd_pointer,entry_2,'predict_Model-'+ai_model)
    sources=''
    if  result !=False and use_saved==True: 
        predict_reaction_results=result
        smiles=predict_reaction_results['response']['payload']['attempts'][0]['smiles']
    
        confidence=predict_reaction_results['response']['payload']['attempts'][0]['confidence']
        x_y=smiles.split('>>')[1]
        source=[]
        
        if cmd_pointer.notebook_mode == True:
            from IPython.display import display,Markdown
            
            display(Markdown('***Saved Result*** '))
            
            display(rxn_helper.output('<success>Smiles:</success>    '+smiles,cmd_pointer ))
                    
            
            for x in source:
                if len(sources) > 0:
                    sources=sources+' + '+x
                else: 
                    sources=x
            
            display(rxn_helper.output('<success>Reaction:</success> '+sources+'    ---->    '+x_y ,cmd_pointer))
            display(rxn_helper.output('<success>Confidence:</success> '+str(confidence),cmd_pointer))
                    
            return get_reaction_from_smiles(smiles)
        else:
            print(f'\nSmiles: {smiles}')
            for x in source:
                if len(sources) > 0:
                    sources=sources+' + '+x
                else: 
                    sources=x
            print('Reaction: '+sources+'    ---->    '+x_y )
            print(f'Confidence: {confidence}')
            print(smiles)
            return True
        
   
    try:
        if predict_id==None and ai_model==None:
            
            predict_reaction_response = rxn4chemistry_wrapper.predict_reaction(molecule)
        elif predict_id==None:
            predict_reaction_response = rxn4chemistry_wrapper.predict_reaction(molecule,ai_model=ai_model)
        else:
            predict_reaction_response = rxn4chemistry_wrapper.predict_reaction(molecule,prediction_id=predict_id)
    except BaseException as e:
        print(e)
        raise('\n'+"Prediction call raised an error: "+str(e)+'\n')
    status=False
    retries=0
    while status==False:
            try:           
                if predict_id==None and ai_model==None:
                    predict_reaction_response = rxn4chemistry_wrapper.predict_reaction(molecule)
                elif predict_id==None:
                    predict_reaction_response = rxn4chemistry_wrapper.predict_reaction(molecule,ai_model=ai_model)
                else:
                    predict_reaction_response = rxn4chemistry_wrapper.predict_reaction(molecule,prediction_id=predict_id)
                status=True
            
            except BaseException as e: 
                retries=retries+1
                sleep(2)
                if retries > 4:
                    raise BaseException("Server unresponsive"+str(e))
    sleep(2)
    predict_reaction_results={}
    retries=0
    while 'response' not in predict_reaction_results:
           
            try:
                predict_reaction_results = rxn4chemistry_wrapper.get_predict_reaction_results(predict_reaction_response['prediction_id'])
                if 'response' not in predict_reaction_results:
                    sleep(2)
            except BaseException as e: 
                retries=retries+1
                if retries > 4:
                    raise BaseException("Server unresponsive"+str(e))


    
    smiles=predict_reaction_results['response']['payload']['attempts'][0]['smiles']
    
    confidence=predict_reaction_results['response']['payload']['attempts'][0]['confidence']
    rxn_helper.save_to_results_cache(cmd_pointer,smiles.split('>>')[0].split('.'),predict_reaction_results,'predict_Model-'+ai_model)
            
    source=[]
    for i in predict_reaction_results['response']['payload']['attempts'][0]['smiles'].split(">>")[0].split('.'):
        source.append(i)
    x_y=predict_reaction_results['response']['payload']['attempts'][0]['smiles'].split('>>')[1]
        #print(helper._smiles_to_iupac(reaction_prediction["smiles"].split(".")[0]))
   
    if cmd_pointer.notebook_mode == True:
        from IPython.display import display,Markdown
        
        display(Markdown('***Generated Result*** '))
         
        display(rxn_helper.output('<success>Smiles:</success>    '+smiles,cmd_pointer ))
                
        sources=''
        for x in source:
            if len(sources) > 0:
                sources=sources+' + '+x
            else: 
                sources=x
        
        display(rxn_helper.output('<success>Reaction:</success> '+sources+'    ---->    '+x_y ,cmd_pointer))
        display(rxn_helper.output('<success>Confidence:</success> '+str(confidence),cmd_pointer))
                
        return get_reaction_from_smiles(predict_reaction_results['response']['payload']['attempts'][0]['smiles'])
    else:
        print(f'\nSmiles: {smiles}')
        for x in source:
            if len(sources) > 0:
                sources=sources+' + '+x
            else: 
                sources=x
        print('Reaction: '+sources+'    ---->    '+x_y )
        print(f'Confidence: {confidence}')
        print(smiles)
        return True
       


def confirm_prompt(question: str) -> bool:
    import readline
    reply = None
    while reply not in ("y", "n"):
        reply = input(f"{question} (y/n): ").casefold()
        readline.remove_history_item(readline.get_current_history_length() - 1)
    return (reply == "y")
