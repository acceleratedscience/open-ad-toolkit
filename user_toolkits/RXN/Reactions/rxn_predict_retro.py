

_tableformat = 'simple'

from rxn4chemistry import RXN4ChemistryWrapper
list_of_reactions =[]

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
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from typing import Dict, List
def get_reaction_from_smiles(reaction_smiles: str) -> Chem.rdChemReactions.ChemicalReaction:
    return AllChem.ReactionFromSmarts(reaction_smiles, useSmiles=True)
# https://cps.foc-deepsearch.zurich.ibm.com/projects/1234567890abcdefghijklmnopqrstvwyz123456/library

def collect_reactions_from_retrosynthesis(tree: Dict)  -> List[str] :

    reactions = []
    global list_of_reactions
    
    if 'children' in tree and len(tree['children']):
        reactions.append(
            AllChem.ReactionFromSmarts('{}>>{}'.format(
                '.'.join([node['smiles'] for node in tree['children']]),
                tree['smiles']
            ), useSmiles=True)
        )
    
    for node in tree['children']:
        reactions.extend(collect_reactions_from_retrosynthesis(node))
    
    return reactions

def collect_reactions_from_retrosynthesis_text(tree: Dict) ->List[str]  :
    
   
    reactions = []
    global list_of_reactions
    
    if 'children' in tree and len(tree['children']):
        reactions.append(
            '{} --->> {}'.format(
                ' + '.join([node['smiles'] for node in tree['children']]),
                tree['smiles']
            ))
    
    for node in tree['children']:
        reactions.extend(collect_reactions_from_retrosynthesis_text(node))
    
    return reactions



def predict_retro(inputs: dict, toolkit_dir, cmd_pointer):
    rxn_helper=get_include_lib(cmd_pointer)
    
    if cmd_pointer.notebook_mode == True:
        import IPython
        from halo import HaloNotebook as Halo
    else:
        from halo import Halo

    
    class Spinner(Halo):
        def __init__(self):
            # Alternative spinners:
            # simpleDotsScrolling, interval=100
            super().__init__(spinner='dots', color='white')

    
    val='val'
    
    availability_pricing_threshold  =   0
    available_smiles                =   None
    exclude_smiles                  =   None
    exclude_substructures           =   None
    exclude_target_molecule         =   False
    fap                             =   0.6
    max_steps                       =   3
    nbeams                          =   10
    pruning_steps                   =   2
    ai_model                        = '2020-07-01'

    product_smiles = inputs["molecule"]
    
    #######################
    if cmd_pointer.notebook_mode == True:
        from IPython.display import display,Markdown
    if not rxn_helper.valid_smiles(str(product_smiles)):
        if cmd_pointer.notebook_mode == True:
            display(Markdown("***Error:*** Invalid Smiles Supplied."))
        else:
            print("Error: Invalid Smiles Supplied.")
        return False
    

    if cmd_pointer.notebook_mode == True:
        import py3Dmol
        style='stick'
        mol = Chem.MolFromSmiles(product_smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
        mblock = Chem.MolToMolBlock(mol)

        view = py3Dmol.view(width=400, height=300)
        view.addModel(mblock, 'mol')
        view.setStyle({style:{}})
        view.zoomTo()
        view.show()
        
        display(Markdown("***Target Molecule:*** "+product_smiles))

    #######################
    
    if "availability_pricing_threshold" in inputs:
        availability_pricing_threshold =int(inputs[ "availability_pricing_threshold"][val])
    
    if "available_smiles"  in inputs: 
        available_smiles = inputs[ "available_smiles"][val]
    
    if "exclude_smiles" in inputs:
        exclude_smiles = inputs[ "exclude_smiles"][val]
    
    if "exclude_substructures" in inputs:
        exclude_substructures = inputs["exclude_substructures"][val]
    
    if "exclude_target_molecule" in inputs:
        if  inputs["exclude_substructures"][val].upper() == "TRUE":
            exclude_target_molecule = True
    
    if "fap" in inputs:
        fap = float(inputs["fap"][val])
    if "max_steps" in inputs:
        max_steps= int(inputs["max_steps"][val])
    if "nbeams" in inputs:
        nbeams = int(inputs["nbeams"][val])
    if "pruning_steps" in inputs:
        pruning_steps = int(inputs["pruning_steps"][val])
    
    if "ai_model" in inputs:
        ai_model= inputs['ai_model'][val]
    
    
   
    rxn4chemistry_wrapper = cmd_pointer.login_settings['client'][cmd_pointer.login_settings['toolkits'].index('RXN') ]
    # Prepare the data query
    from time import sleep
    predict_retro_responses=None
    try:
        print('\n')
        newspin =Spinner()
        
        newspin.start("Starting Retrosynthesis")
        retries=0
        status=False
        while retries <10 and status==False:
            try:
                newspin.text=("Submitting Retrosynthesis ")
                predict_retro_response = rxn4chemistry_wrapper.predict_automatic_retrosynthesis(product_smiles, availability_pricing_threshold=availability_pricing_threshold, available_smiles=available_smiles, exclude_smiles=exclude_smiles, exclude_substructures=exclude_substructures, exclude_target_molecule=exclude_target_molecule, fap=fap, max_steps=max_steps, nbeams=nbeams, pruning_steps=pruning_steps, ai_model=ai_model)
                
                status=True
                
            except Exception as e: 
                sleep(2)
                print(str(e))
                retries=retries+1
                if retries >=10:
                    raise BaseException("Server unresponsive: Unable to submit for processing after 10 retires"+str(e))
     
        if predict_retro_response == None:
            raise BaseException("Server unresponsive: Unable to submit for processing after 10 retires")
           
        retries=0
        if predict_retro_response['response']['payload']['errorMessage'] !=None:

            return predict_retro_response['response']['payload']['errorMessage']
        
        
        status='NEW'
        sleep_wait=5
        while status !="SUCCESS":
            try:
                newspin.text=("Processing Retrosynthesis :"+status)
                predict_automatic_retrosynthesis_results = rxn4chemistry_wrapper.get_predict_automatic_retrosynthesis_results(predict_retro_response['prediction_id'])
                
                if predict_retro_response['response']['payload']==None  :
                    
                    if cmd_pointer.notebook_mode == True:
                        Markdown('\n ***No Result:***  Unable to find path to '+product_smiles)
                        return
                    else: 
                        return  '\n ***No Result:***  Unable to find path to '+product_smiles
                         
                
                
                status=predict_automatic_retrosynthesis_results['status']
                
                sleep(5)
            except BaseException as e: 
                retries=retries+1
                sleep(15)
                newspin.text=("Processing Retrosynthesis :Waiting" )
                if retries > 20:
                    raise BaseException("Server unresponsive: Unable to complete processing for prediction id:'"+predict_retro_response['prediction_id']+"'after 20 retires"+str(e))
    except  BaseException as e:
            newspin.fail('Unable to Process')
            newspin.stop()
            raise BaseException("Unable to complete processing "+str(e)    )#print(predict_automatic_retrosynthesis_results)
            
    reactions_text=[]
    
    try:
        for index, tree in enumerate(predict_automatic_retrosynthesis_results['retrosynthetic_paths']):
            # display(Markdown('Showing path {} with confidence {}:'.format(index, tree['confidence'])))
            for reaction in collect_reactions_from_retrosynthesis_text(tree):
                #display(str(reaction))
                reactions_text.append(str(reaction)  )
    except BaseException as e:
        newspin.fail('Unable to Process')
        newspin.stop()
        raise BaseException ("The following Error message was received while trying to process results:"+str(e))
    i=0
    try:
        newspin.succeed('Finsihed Processing')
        newspin.start()
        newspin.stop()
        results_list=[]
        for index, tree in enumerate(predict_automatic_retrosynthesis_results['retrosynthetic_paths']):
            if cmd_pointer.notebook_mode == True:
                display(Markdown('***Showing path*** {} ***with confidence*** {}:'.format(index, tree['confidence'])))
            else:
                results_list.append('------------------------------------------------------------------------------------------------------------')
                results_list.append('\n \n Showing path {} with confidence {}:'.format(index, tree['confidence']))
            for reaction in collect_reactions_from_retrosynthesis(tree):
                if cmd_pointer.notebook_mode == True:
                    display(Markdown('\n ***Reaction:***  '+reactions_text[i]))
                    
                    display(Chem.Draw.ReactionToImage(reaction))
                else:
                    results_list.append('\n Reaction:  '+reactions_text[i])
                i=i+1
       
        
    except BaseException as e:
        newspin.fail('Unable to Display')
        newspin.stop()
        raise BaseException ("The following Error message was received while trying to display results:"+str(e))
    
    import pandas as pd
    df = pd.DataFrame(results_list)
    from tabulate import tabulate
    return '\n'+tabulate(df, tablefmt=_tableformat,  showindex=False)+'\n'
    