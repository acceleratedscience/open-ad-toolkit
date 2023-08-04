

_tableformat = 'simple'

from rxn4chemistry import RXN4ChemistryWrapper
from rdkit import Chem
from rdkit.Chem import AllChem


def get_reaction_from_smiles(reaction_smiles: str) -> Chem.rdChemReactions.ChemicalReaction:
    return AllChem.ReactionFromSmarts(reaction_smiles, useSmiles=True)

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

# https://cps.foc-deepsearch.zurich.ibm.com/projects/1234567890abcdefghijklmnopqrstvwyz123456/library


def predict_reaction_batch(inputs: dict, toolkit_dir, cmd_pointer):
    
   
    if cmd_pointer.notebook_mode == True:
        import IPython
        from halo import HaloNotebook as Halo
    else:
        from halo import Halo
    
    rxn_helper=get_include_lib(cmd_pointer)
    rxn_helper.sync_up_workspace_name(cmd_pointer)

    name,id = rxn_helper.get_current_project(cmd_pointer)
    if cmd_pointer.notebook_mode == True:
        from IPython.display import display,Markdown
    if name == None and cmd_pointer.api_mode==False:
        if cmd_pointer.notebook_mode==True:
            
            display(Markdown(" No current RXN project selected ,`set rxn project <project name>` to set your project before proceeding."))
            display(Markdown("Select from the Below Projects or create a new."))
            display(rxn_helper.get_all_projects(cmd_pointer)[['name','description']])
        else:
            print(" No current RXN project selected ,`set rxn project <project name>` to set your project before proceeding. ")
            print(" Select from the Below Projects or create a new.")
            print(rxn_helper.get_all_projects(cmd_pointer)[['name','description']])
        return False
    
 
    class Spinner(Halo):
        def __init__(self):
            # Alternative spinners:
            # simpleDotsScrolling, interval=100
            super().__init__(spinner='dots', color='white')
   
    
    if 'from_list' not in inputs['from_source'][0]:
        print('lists are the only source type implemented current, dataframe and file are soon to be implemented')
        return True
    else:
        from_list= inputs['from_source'][0]['from_list']
    
   
    
    
    newspin =Spinner()
  
    newspin.start("Starting Prediction")
   
    
    
    val='val'
    if "ai_model" in inputs:
        ai_model= inputs['ai_model'][val]
    else:
        ai_model="2020-08-10"
    new_from_list=[]
    cached_results=[]
    if 'use_saved'in inputs:
         use_saved = True
    else:
        use_saved=False
    from rdkit import Chem

    for entry in from_list:
        error_list=[]
        for i in entry.split('.'):
            if not rxn_helper.valid_smiles(str(i)):
                error_list.append(i)
            
        if len(error_list)>0:
            import pandas as pd
            df = pd.DataFrame(error_list,columns=['smiles'])
            if cmd_pointer.notebook_mode == True:
                
                display(Markdown("***Error:*** The following invalid were Smiles Supplied:"))
                display(df)
                display(Markdown("***info:*** This reaction will be skipped  "+entry+" "))
            else:
                print("Error: The following invalid were Smiles Supplied:")
                
                print(df)
                print("Info: This reaction will be skipped "+entry+"/n/n")
            
            continue         
       
        entry_2=[]
        dot=''
        for i in entry.split('.'):
            entry_2.append(Chem.MolToSmiles(Chem.MolFromSmiles(i),canonical=True))
            dot='.'
       
        result = rxn_helper.retrieve_cache(cmd_pointer,entry_2,'predict_batch_Model-'+ai_model)
        
        
        
        if  result !=False and use_saved==True:
            cached_results.append(result)
        else:
            new_from_list.append(entry)    

    results_list=[]

    for reaction_prediction in cached_results:
        
        source=[]
        for i in reaction_prediction['smiles'].split(">>")[0].split('.'):

            source.append(i)
        x_y=reaction_prediction['smiles'].split('>>')[1]
        
        if cmd_pointer.notebook_mode == True:
            from IPython.display import display,Markdown
            display(Markdown('***Saved Result*** '))
          
            display(rxn_helper.output(f'<success>Smiles:</success> {reaction_prediction["smiles"]}',cmd_pointer))
            sources=''
            for x in source:
                if len(sources) > 0:
                    sources=sources+' + '+x
                else: 
                    sources=x
            
            display(rxn_helper.output(f'<success>Reaction:</success> '+sources+'    ---->    '+x_y ,cmd_pointer))
            display(rxn_helper.output(f'<success>Confidence:</success> {reaction_prediction["confidence"]}',cmd_pointer))
        
            display(get_reaction_from_smiles(reaction_prediction['smiles']))
        else:
            results_list.append('CACHED RESULT')
            results_list.append(f'Smiles: {reaction_prediction["smiles"]}')
            sources=''
            for x in source:
                if len(sources) > 0:
                    sources=sources+' + '+x
                else: 
                    sources=x
            results_list.append('Reaction: '+sources+'    ---->    '+x_y )
            results_list.append(f'Confidence: {reaction_prediction["confidence"]}')
            results_list.append('____________________________________________________')
    from time import sleep
    if len(new_from_list) > 0:
        from_list=new_from_list
        rxn4chemistry_wrapper = cmd_pointer.login_settings['client'][cmd_pointer.login_settings['toolkits'].index('RXN') ]
        retries=0
        status=False
        while status==False:
            
           
            try:
                newspin.text=("Processing Prediction" )
                
                predict_reaction_batch_response = rxn4chemistry_wrapper.predict_reaction_batch(from_list)
                sleep(2)
                status=True
            except BaseException as e: 
                retries=retries+1
                if retries > 4:
                    newspin.fail('Unable to Process')
                    newspin.stop()
                    raise BaseException("Server unresponsive"+str(e))
        
       
       
        
        
        
        retries=0
        
        #reaction_predictions=rxn4chemistry_wrapper.get_predict_reaction_batch_results(predict_reaction_batch_response['task_id'])
        reaction_predictions={}
        while 'predictions' not in reaction_predictions:
           
            try:
                newspin.text=("Processing Prediction" )
                
                reaction_predictions=rxn4chemistry_wrapper.get_predict_reaction_batch_results(predict_reaction_batch_response['task_id'])
                if 'predictions' not in reaction_predictions:
                    sleep(3)
            except BaseException as e: 
                retries=retries+1
                if retries > 10:
                    newspin.fail('Unable to Process')
                    newspin.stop()
                    raise BaseException("Server unresponsive"+str(e))
        
          
        for reaction_prediction in reaction_predictions['predictions']:
          
            rxn_helper.save_to_results_cache(cmd_pointer,reaction_prediction['smiles'].split('>>')[0].split('.'),reaction_prediction,'predict_batch_Model-'+ai_model)
            source=[]
            for i in reaction_prediction['smiles'].split(">>")[0].split('.'):
                source.append(i)
            x_y=reaction_prediction['smiles'].split('>>')[1]
            
            if cmd_pointer.notebook_mode == True:
                from IPython.display import display,Markdown
                display(Markdown('***Generated Result*** '))
                display(rxn_helper.output(f'<success>Smiles:</success> {reaction_prediction["smiles"]}',cmd_pointer))
                sources=''
                for x in source:
                    if len(sources) > 0:
                        sources=sources+' + '+x
                    else: 
                        sources=x
                display(rxn_helper.output('<success>Reaction:</success> '+sources+'    ---->    '+x_y ,cmd_pointer))
                display(rxn_helper.output(f'<success>Confidence:</success> {reaction_prediction["confidence"]}',cmd_pointer))
                display(get_reaction_from_smiles(reaction_prediction['smiles']))
            else:
                results_list.append(f'Smiles: {reaction_prediction["smiles"]}')
                sources=''
                for x in source:
                    if len(sources) > 0:
                        sources=sources+' + '+x
                    else: 
                        sources=x
                results_list.append('Reaction: '+sources+'    ---->    '+x_y )
                results_list.append(f'Confidence: {reaction_prediction["confidence"]}')
                results_list.append('____________________________________________________')
    newspin.succeed('Finsihed Processing')
    newspin.start()
    newspin.stop()        
    if cmd_pointer.notebook_mode != True:
        import pandas as pd
        df = pd.DataFrame(results_list)
        from tabulate import tabulate
        return '\n'+tabulate(df, tablefmt=_tableformat,  showindex=False)+'\n'
    
    