

_tableformat = 'simple'

from rxn4chemistry import RXN4ChemistryWrapper
from rdkit import Chem
from rdkit.Chem import AllChem

from time import sleep

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


def predict_reaction_batch_topn(inputs: dict, toolkit_dir, cmd_pointer):
    top_n=5
    rxn_helper=get_include_lib(cmd_pointer)
    rxn_helper.sync_up_workspace_name(cmd_pointer)
    name,id= rxn_helper.get_current_project(cmd_pointer)

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
    
    class Spinner(Halo):
        def __init__(self):
            # Alternative spinners:
            # simpleDotsScrolling, interval=100
            super().__init__(spinner='dots', color='white')
    
    if cmd_pointer.notebook_mode == True:
                from IPython.display import display,Markdown
   

    val='val'
    if 'topn' in inputs:
        
        if int(inputs['topn'][val]) > 5:
            top_n=int(inputs['topn'][val])
    
    if "ai_model" in inputs:
        ai_model= inputs['ai_model'][val]
    else:
        ai_model="2020-08-10"

    rxn_helper = get_include_lib(cmd_pointer)
    ###################################################################################################
   # getting our input source for the reactions
    
    if 'from_list'  in inputs['from_source'][0]:
        from_list= inputs['from_source'][0]['from_list']
    elif 'from_dataframe' in inputs:
        try:
            react_frame = cmd_pointer.api_variables[inputs['from_dataframe']]
            from_list=rxn_helper.get_column_as_list_from_dataframe(react_frame,'reactions')
            if from_list == []:
                raise BaseException
        except BaseException as err:
           print ("Could not load valid list from dataframe column 'reactions' ")
           return True
    elif 'from_file' in inputs:
        from_file= inputs['from_file']
        try:
            react_frame = rxn_helper.get_dataframe_from_file(cmd_pointer,from_file)
        
            from_list=rxn_helper.get_column_as_list_from_dataframe(react_frame,'reactions')
            if from_list == []:
                raise BaseException
        except BaseException as err:
           
           print ("Could not load valid list from file column 'reactions' ")
           return True
    
    #print(inputs['use_saved'])
    if 'use_saved'in inputs:
         use_saved = True
    else:
        use_saved=False
    new_from_list=[]
    cached_results=[]
    
    from rdkit import Chem
    new_cannonical_list=[]
    cached_cannonical_list=[]
    results_list=[]
    
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
                display(Markdown("***info:*** This rection will be skipped: "+entry+" "))
            else:
                print("Error: The following invalid were Smiles Supplied:")

                print(df)
                print("Info: This reaction will be skipped: "+entry+"/n/n")
            
            continue    
        entry_2=[]
        dot=''
   
        for i in entry.split('.'):
            entry_2.append(Chem.MolToSmiles(Chem.MolFromSmiles(i),canonical=True))
            dot='.'
        
        result = rxn_helper.retrieve_cache(cmd_pointer,entry_2,'predict_batch_topn'+str(top_n)+"_model"+ai_model)
        if  result !=False and use_saved==True:
            
            cached_results.append(result)
            cached_cannonical_list.append(entry_2)
        else:
            new_from_list.append(entry.split('.'))    
            new_cannonical_list.append(entry_2)
    reaction_no=0
    
    for reaction_predictions in cached_results:
        if cmd_pointer.notebook_mode == True:
                 display(Markdown(f"\n\n***Saved Reaction results for:*** {cached_cannonical_list[reaction_no]}"))
        else:
            results_list.append(f"\n\nSaved Reaction results for: {cached_cannonical_list[reaction_no]}")
        reaction_no=reaction_no+1
        for j, prediction in enumerate(reaction_predictions["results"], 1):
                product_smiles = ".".join(prediction["smiles"])
                confidence = prediction["confidence"]
                if cmd_pointer.notebook_mode == True:
                    display(rxn_helper.output(cmd_pointer,f'<success>         Product(s) {j} {product_smiles}, With confidence {confidence}</success>'),cmd_pointer)
                    #display(Markdown(f'         Product(s) {j} {product_smiles}, With confidence {confidence}'))
                else:          
                    results_list.append(f'        Product(s) {j}: {product_smiles}, with confidence {confidence}')
        results_list.append('____________________________________________________')    
    if len(new_from_list) > 0:
        
    
        val='val'
        
        newspin =Spinner()
        newspin.start("Starting Prediction")
        from_list=new_from_list
       
        retries=0
        status=False
        rxn4chemistry_wrapper = cmd_pointer.login_settings['client'][cmd_pointer.login_settings['toolkits'].index('RXN') ]
        while status==False:
            try:
                newspin.text=("Processing Prediction" )
                sleep(2)
                
                predict_rection_batch_response = rxn4chemistry_wrapper.predict_reaction_batch_topn(precursors_lists=new_from_list,topn=top_n,ai_model=ai_model,)
                status=True
            except BaseException as e: 
                retries=retries+1
                if retries > 4:
                    newspin.fail('Unable to Process')
                    newspin.stop()
                    raise BaseException("Server unresponsive"+str(e))
        
        x={}
        retries=0
        while 'predictions' not in x:
           
            try:
                newspin.text=("Processing Prediction" )
               
                x=rxn4chemistry_wrapper.get_predict_reaction_batch_topn_results(predict_rection_batch_response['task_id'])
                if 'predictions' not in x:
                    sleep(3)
            except BaseException as e: 
                retries=retries+1
                if retries > 10:
                    newspin.fail('Unable to Process')
                    newspin.stop()
                    raise BaseException("Server unresponsive"+str(e))
        

                
        
    

        reaction_no=0
        for i, reaction_predictions in enumerate(x['predictions'], 1):
            if cmd_pointer.notebook_mode == True:
                 display(Markdown(f'\n\n***Outcomes for reaction:***   {new_cannonical_list[reaction_no]}:'))
            else:          
                results_list.append(f'\n\n Outcomes for reaction   {new_cannonical_list[reaction_no]}:')
            rxn_helper.save_to_results_cache(cmd_pointer,new_cannonical_list[reaction_no],reaction_predictions,'predict_batch_topn'+str(top_n)+"_model"+ai_model)
            reaction_no=reaction_no+1
            for j, prediction in enumerate(reaction_predictions["results"], 1):
                product_smiles = ".".join(prediction["smiles"])
                confidence = prediction["confidence"]
                if cmd_pointer.notebook_mode == True:
                    display(rxn_helper.output(f'<success>         Product(s) {j} {product_smiles}, With confidence {confidence}</success>',cmd_pointer))
                else:          
                    results_list.append(f'        Product(s) {j}: {product_smiles}, with confidence {confidence}')
            results_list.append('____________________________________________________')

        newspin.succeed('Finished Processing')
    newspin.start()
    newspin.stop()        
    if cmd_pointer.notebook_mode != True:
        import pandas as pd
        df = pd.DataFrame(results_list)
        from tabulate import tabulate
        return '\n'+tabulate(df, tablefmt=_tableformat,  showindex=False)+'\n'
    
        
