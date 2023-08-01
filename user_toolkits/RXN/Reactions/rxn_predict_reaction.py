

_tableformat = 'simple'

from rxn4chemistry import RXN4ChemistryWrapper
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
    
    rxn4chemistry_wrapper = cmd_pointer.login_settings['client'][cmd_pointer.login_settings['toolkits'].index('RXN') ]
    # Prepare the data query
    rxn_helper = get_include_lib(cmd_pointer) 
    molecule = inputs["molecule"]
    val='val'
    predict_id=None
    ai_model=None
    if 'prediction_id' in inputs:
       predict_id = inputs["prediction_id"][val]
    if 'ai_model' in inputs:
        ai_model = inputs["ai_model"][val]
    if cmd_pointer.notebook_mode == True:
        from IPython.display import display,Markdown
    error_list=[]
    
    
    
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
        print(f'Smiles: {smiles}')
        for x in source:
            if len(sources) > 0:
                sources=sources+' + '+x
            else: 
                sources=x
        print('Reaction: '+sources+'    ---->    '+x_y )
        print(f'Confidence: {confidence}')
        print(smiles)
        return True

    




    import pandas as pd
    df = pd.DataFrame(x)
    #df['createdOn'] = pd.to_datetime(df['createdOn'], format='mixed')
    #df['modifiedOn'] = pd.to_datetime(df['modifiedOn'], format='mixed')
    #df= df[['name','description','attempts']]
    
    if cmd_pointer.notebook_mode == True:
        return df

    from tabulate import tabulate
    print('\n'+tabulate(df, tablefmt=_tableformat, headers="keys", showindex=False)+'\n')
    
    return True
    if result is None:
        print("Search returned no result")
        return None

    if cmd_pointer.notebook_mode == True:
        import pandas as pd
        df = pd.DataFrame(results_table)
        import numpy as np
        if 'save_as' in inputs:
            df.to_csv(cmd_pointer.workspace_path(cmd_pointer.settings['workspace'].upper()) + "/" + inputs['results_file'] + '.csv',)
        return (df. replace(np.nan, '', regex=True))

    else:
        from tabulate import tabulate
        import pandas as pd
        import numpy as np
        collectives = pd.DataFrame(results_table)

        if 'save_as' in inputs:

            collectives.to_csv(
                cmd_pointer.workspace_path(
                    cmd_pointer.settings['workspace'].upper()) +
                "/" +
                inputs['results_file'] +
                '.csv',
                index=False)

       


def confirm_prompt(question: str) -> bool:
    import readline
    reply = None
    while reply not in ("y", "n"):
        reply = input(f"{question} (y/n): ").casefold()
        readline.remove_history_item(readline.get_current_history_length() - 1)
    return (reply == "y")
