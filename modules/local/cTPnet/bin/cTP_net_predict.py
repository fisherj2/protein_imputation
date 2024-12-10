import numpy as np
import torch.nn as nn
import pandas as pd
import torch.nn.functional as F
import torch
import sys

import inspect
import os
script_directory = os.path.dirname(os.path.abspath(
  inspect.getfile(inspect.currentframe())))
 
print(script_directory)


import importlib.util
 
# specify the module that needs to be 
# imported relative to the path of the 
# module
spec=importlib.util.spec_from_file_location("module", script_directory + "/cTP_net_network.py")
 
# creates a new module based on spec
cTP_net_network=importlib.util.module_from_spec(spec)
 
# executes the module in its own namespace
# when a module is imported or reloaded.
spec.loader.exec_module(cTP_net_network)

# from . import cTP_net_network


def modified_predict(X,model_file_path,d=12,protein_names='',rna_dim=0):
    """
    Main function: imputation of 12 or 24, (or custom) surface protein abundance from single
    cell RNA-seq matrix.

    Parameters
    ----------
    arg1 : pandas DataFrame
        X: A pandas data.frame with dimension DxN, where D equals to 12611 as 
        the number of genes we trained the neural network; and N is the number 
        of cells. N>=1
    arg2 : str
        path to the trained cTPnet pytorch model for the prediction
    arg3 : str
        vector of protein names used in training data
    arg3 : int
        number of RNA features in training data

    Returns
    -------
    pandas DataFrame
        Impuated surface protein abundance DataFrame, with rows for proteins and
        columns for cells

    """
    


    print('custom class defined')
    print('python package loaded')
    sys.stdout.flush()
    X=X.apply(lambda x: np.log((x*10000.0/sum(x))+1.0), axis=0)
    print('convert input')
    sys.stdout.flush()
    cells=X.columns
    if d==12:
        proteins=['CD45RA', 'CD19', 'CD14', 'CD8', 'CD16', 'CD3', 'CD11c', 'CD4','CD56','CD34','CD2','CD57']
        net=cTP_net_network.Net12()
    elif d==24:
        proteins=['CD3', 'CD4', 'CD45RA',  'CD16', 'CD14','CD11c', 'CD19','CD8','CD34', 'CD56','CD57','CD2','CD11a','CD123','CD127-IL7Ra','CD161','CD27','CD278-ICOS','CD28','CD38','CD45RO','CD69','CD79b','HLA.DR']
        net=cTP_net_network.Net24()
    elif d=='custom':
      
        #build required Net class
        class Net(nn.Module):
            def __init__(self):
                super(Net, self).__init__()
                #self.fc1 = nn.Linear(14505, 1000)
                self.fc1 = nn.Linear(rna_dim, 1000)
                self.fc2 = nn.Linear(1000, 256)
                self.fc3 = nn.ModuleDict({})
                for p in protein_names:
                	self.fc3[p]=nn.Linear(256, 64)
                self.fc4 = nn.ModuleDict({})
                for p in protein_names:
                	self.fc4[p]=nn.Linear(64, 1)
            def forward(self, x):
                x = F.relu(self.fc1(x))
                x = F.relu(self.fc2(x))
                outputs={}
                for p in protein_names:
                	outputs[p]=self.fc4[p](F.relu(self.fc3[p](x)))
                return outputs
              
        proteins=protein_names
        net=Net()
    else:
        raise Exception('Protein output dimension has to be 12 or 24, or be set to use custom Net class definition')
    
    print('init network')
    sys.stdout.flush()
    net.load_state_dict(torch.load(model_file_path))
    print('load network')
    sys.stdout.flush()
    X=torch.t(torch.tensor(X.values))
    X=X.type(torch.FloatTensor)
    
    if d=='custom':
      pred_dict=net(X)
      y_pred=list(pred_dict.values())
    else:
      y_pred = list(net(X))
    
    
    print('imputation')
    sys.stdout.flush()
    y_pred=torch.transpose(torch.stack(y_pred),0,1).view(X.shape[0],-1)
    y_pred=pd.DataFrame(y_pred.detach().numpy().T)
    y_pred.index=proteins
    y_pred.columns=cells
    print('convert format')
    sys.stdout.flush()
    print('python done')
    return(y_pred)

