
import scanpy as sc
import scipy.sparse as sp
import numpy as np
import torch
import argparse
import anndata as ad
import pandas as pd
from anndata import AnnData, read_h5ad
from scanpy import read
from torch.utils.data import DataLoader, TensorDataset
from torch import tensor
from torch.cuda import is_available
from sciPENN.sciPENN_API import sciPENN_API
import os 

#set seed
seed = 123
torch.manual_seed(seed)
torch.cuda.manual_seed(seed) 
torch.cuda.manual_seed_all(seed)  # if you are using multi-GPU
np.random.seed(seed)

#setup arg parse
parser = argparse.ArgumentParser( prog = 'Script to train scipenn on scRNAseq data')


parser.add_argument('-d', '--basedir', required=True, help="pipeline base directory")
parser.add_argument('-l', '--launchdir', required=True, help="pipeline launch directory")
parser.add_argument('-b','--bench',  help='<Required> Set flag for benchmarking', required=True)
parser.add_argument('-f','--files', nargs='+', help='<Required> Set input files', required=True)


args = parser.parse_args()
dobenchmark = args.bench
input_files = args.files
print(input_files)


#if length is one, probably need to split
if len(input_files)==1 or isinstance(input_files, (str)):
    input_files=input_files[0]
    input_files=input_files.split('.csv')
    input_files=[s + '.csv' for s in input_files]

#check output dir exists
if not os.path.exists(args.launchdir + "/output/sciPENN"):
    os.makedirs(args.launchdir + "/output/sciPENN")

#need to clean inputs of wrong character
characters_to_remove = ['[', ']', ',']
translation_table = str.maketrans('', '', ''.join(characters_to_remove))



metadata_file_train =  [file for file in input_files if 'metadata_train' in file][0]
metadata_file_train = metadata_file_train.translate(translation_table)
metadata_file_test =  [file for file in input_files if 'metadata_test' in file][0]
metadata_file_test = metadata_file_test.translate(translation_table)

rna_train_data_file = [file for file in input_files if 'training_data_rna_raw' in file][0]
rna_train_data_file = rna_train_data_file.translate(translation_table)

prot_train_data_file = [file for file in input_files if 'training_data_prot_raw' in file][0]
prot_train_data_file = prot_train_data_file.translate(translation_table)

rna_test_data_file =  [file for file in input_files if 'testing_data_rna_raw' in file][0]
rna_test_data_file = rna_test_data_file.translate(translation_table)

#remove whitespace
metadata_file_train=metadata_file_train.strip()
metadata_file_test=metadata_file_test.strip()
rna_train_data_file=rna_train_data_file.strip()
prot_train_data_file=prot_train_data_file.strip()
rna_test_data_file=rna_test_data_file.strip()

#metadata_file =  '/scratch/jfisher2/protein_prediction/test_output/training_files/metadata.csv'
#rna_train_data_file = '/scratch/jfisher2/protein_prediction/test_output/training_files/training_data_rna_raw.csv'
#prot_train_data_file = '/scratch/jfisher2/protein_prediction/test_output/training_files/training_data_prot_raw.csv'
#rna_test_data_file =  '/scratch/jfisher2/protein_prediction/test_output/testing_files/testing_data_rna_raw.csv'

adata_rna_train = ad.read_csv(rna_train_data_file, delimiter=',', first_column_names = True).transpose()
adata_prot_train = ad.read_csv(prot_train_data_file, delimiter=',', first_column_names = True).transpose()
adata_rna_test = ad.read_csv(rna_test_data_file, delimiter=',', first_column_names = True).transpose()

# # # read in as pandas
meta_train = pd.read_csv(metadata_file_train, index_col=0)
meta_test = pd.read_csv(metadata_file_test, index_col=0)

#add metadata
ind=meta_train.index.isin(adata_rna_train.obs.index)
adata_rna_train.obs['donor']=meta_train.loc[ind,'donor']
adata_rna_train.obs['celltype']=meta_train.loc[ind,'celltype']

ind=meta_train.index.isin(adata_prot_train.obs.index)
adata_prot_train.obs['donor']=meta_train.loc[ind,'donor']
adata_prot_train.obs['celltype']=meta_train.loc[ind,'celltype']

ind=meta_test.index.isin(adata_rna_test.obs.index)
adata_rna_test.obs['donor']=meta_test.loc[ind,'donor']
adata_rna_test.obs['celltype']=meta_test.loc[ind,'celltype']

#fix metadata, convert to categorical dtype
adata_rna_train.obs['donor'] = adata_rna_train.obs['donor'].astype('category')
adata_prot_train.obs['donor'] = adata_prot_train.obs['donor'].astype('category')
adata_rna_test.obs['donor'] = adata_rna_test.obs['donor'].astype('category')
adata_rna_train.obs['celltype'] = adata_rna_train.obs['celltype'].astype('category')
adata_prot_train.obs['celltype'] = adata_prot_train.obs['celltype'].astype('category')
adata_rna_test.obs['celltype'] = adata_rna_test.obs['celltype'].astype('category')

adata_rna_train.var['features'] = pd.DataFrame(adata_rna_train.var.index, index=adata_rna_train.var.index)
adata_rna_train.var['features'] = adata_rna_train.var['features'].astype('category')
adata_prot_train.var['features'] = pd.DataFrame(adata_prot_train.var.index, index=adata_prot_train.var.index)
adata_prot_train.var['features'] = adata_prot_train.var['features'].astype('category')
adata_rna_test.var['features'] = pd.DataFrame(adata_rna_test.var.index, index=adata_rna_test.var.index)
adata_rna_test.var['features'] = adata_rna_test.var['features'].astype('category')

#convert to sparse matrix
adata_rna_train.X = sp.csc_matrix(adata_rna_train.X)
adata_prot_train.X = sp.csc_matrix(adata_prot_train.X)
adata_rna_test.X = sp.csc_matrix(adata_rna_test.X)



#let scipenn function do normalisation as part of training
sciPENN = sciPENN_API(gene_trainsets = [adata_rna_train], protein_trainsets = [adata_prot_train],
                      gene_test = adata_rna_test, train_batchkeys = ['donor'], test_batchkey = 'donor',
                      type_key = 'celltype', gene_list = [], select_hvg = True, cell_normalize = True, log_normalize = True, gene_normalize = True, min_cells = 1, min_genes = 1, 
                      batch_size = 128, val_split = 0.1, use_gpu = True)
                      

sciPENN.train(weights_dir = "output/sciPENN/scipenn_weights", quantiles = [0.1, 0.25, 0.75, 0.9], n_epochs = 10000, ES_max = 12, decay_max = 6, decay_step = 0.1, lr = 10**(-3),load = True)
 

predicted_test = sciPENN.predict()


    
predicted_test.write_h5ad(
    args.launchdir + "/output/sciPENN/sciPENN_prediction.h5ad"
)

if dobenchmark=='true':
    #get benchmarking prefix

    #make sure looking at filename only, no extended file path
    basename=os.path.basename(rna_train_data_file)
    prefix= basename.replace("_training_data_rna_norm.csv", "") 
    print(prefix)
    np.savetxt(args.launchdir + "/output/sciPENN/" +prefix + "sciPENN_prediction.csv", predicted_test.X, delimiter=",")
    np.savetxt(prefix + "sciPENN_prediction.csv", predicted_test.X, delimiter=",")
else:
    #save for later
    np.savetxt(args.launchdir + "/output/sciPENN/sciPENN_prediction.csv", predicted_test.X, delimiter=",")

    #pass to pipeline
    np.savetxt("sciPENN_prediction.csv", predicted_test.X, delimiter=",")
