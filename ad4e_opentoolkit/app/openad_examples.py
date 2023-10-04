import os,sys
import shutil
from  ad4e_opentoolkit.app.global_var_lib import _repo_dir
destination=None
def openad_create_examps():
    
    try:
        try:
            shutil.rmtree(os.path.expanduser('~/openad_notebooks'),ignore_errors=True)
            #os.mkdir(os.path.expanduser('~/openad_notebooks'))
        except:
            pass
    
        shutil.copytree(_repo_dir+'/../notebooks', os.path.abspath(os.path.expanduser('~/openad_notebooks')),dirs_exist_ok=True,symlinks=True)
        
    except Exception as e:
        print("Error: Unable to copy magic file, destination may not exist")
        print(e)
