from IPython.core.magic import (Magics,magics_class,line_magic,cell_magic,line_cell_magic,needs_local_scope)
from ipywidgets import HTML
from ipywidgets.widgets import Output
from IPython.core.interactiveshell import InteractiveShell
InteractiveShell.ast_node_interactivity = "all"
import os
import time
import sys
sys.path.insert(0, '../')
os.sys.path.append(os.path.dirname(os.path.abspath('./')))
module_path = os.path.abspath(os.path.join('..'))
import tabulate
from halo import HaloNotebook as Halo
import main
handle_cache       =   {'toolkits':[],'toolkits_details':[],'toolkits_api':[],'client':[],'expiry':[],'session_vars':[]}
context_cache      =  {'workspace':None,'toolkit':None}

@magics_class
class AD(Magics):
    
    @needs_local_scope
    @line_cell_magic
    def adccl(self,line,cell=None, local_ns=None):
        import pandas
        import re
        api_variable={}
        #line=line.replace("\n",'')
        #line=line.replace("\t",'')
        #print(line)
        line_list= line.split()
        x=len(line_list)
        i=1
        if x > 1:
            while i < x:
                if line_list[i-1].upper() == 'DATAFRAME':
                    try:
                        
                        df = eval(line_list[i])
                        
                        if isinstance(df, pandas.DataFrame):
                            api_variable[line_list[i]]=df
                    except:
                        pass
                i+=1
     
        return main.api_remote(line,handle_cache,context_cache,api_variable)
      
ip = get_ipython()
ip.register_magics(AD)