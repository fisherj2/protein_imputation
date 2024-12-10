
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import pickle
from time import time
from scipy.stats import spearmanr, gamma, poisson
import anndata as ad
from anndata import AnnData, read_h5ad
import scanpy as sc
from scanpy import read
import scipy.sparse as sp
import torch
from torch.utils.data import DataLoader, TensorDataset
from torch import tensor
from torch.cuda import is_available
from scMMT.scMMT_API import  preprocess
from modified_scMMT_API import scMMT_API
from sklearn.metrics import f1_score, accuracy_score

import argparse
import os

import contextlib
with contextlib.redirect_stdout(None):
    from scMMT.scMMT_API import  preprocess

import os


#setup arg parse
parser = argparse.ArgumentParser( prog = 'Script to train scMMT on scRNAseq data')
parser.add_argument('-d', '--basedir', required=True, help="pipeline base directory")
parser.add_argument('-b','--bench',  help='<Required> Set flag for benchmarking', required=True)
parser.add_argument('-f','--files', nargs='+', help='<Required> Set flag', required=True)

args = parser.parse_args()
dobenchmark = args.bench
input_files = args.files

#if length is one, probably need to split
if len(input_files)==1 or isinstance(input_files, (str)):
    print('adjusting input format')
    if isinstance(input_files, (list)):
        input_files=input_files[0]
    
    input_files=input_files.split('.csv')
    input_files=[s + '.csv' for s in input_files]

print(input_files)

#check output dir exists
if not os.path.exists(args.basedir + "/output/scMMT"):
    os.makedirs(args.basedir + "/output/scMMT")

#need to clean inputs of wrong character
characters_to_remove = ['[', ']', ',']
translation_table = str.maketrans('', '', ''.join(characters_to_remove))


metadata_file_train =  [file for file in input_files if 'metadata_train' in file][0]
metadata_file_train = metadata_file_train.translate(translation_table)
metadata_file_test =  [file for file in input_files if 'metadata_test' in file][0]
metadata_file_test = metadata_file_test.translate(translation_table)

rna_train_data_file = [file for file in input_files if 'training_data_rna_raw' in file][0]
rna_train_data_file = rna_train_data_file.translate(translation_table)

prot_train_data_file = [file for file in input_files if 'training_data_prot_norm' in file][0]
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

adata_rna_test.write_h5ad(
    "/scratch/jfisher2/adata_rna_test.h5ad"
)

scMMT = scMMT_API(    gene_trainsets = [adata_rna_train], protein_trainsets = [adata_prot_train], gene_test = adata_rna_test, 
                      train_batchkeys = ['donor'], test_batchkey = 'donor',
                      log_normalize = True,            # Is scRNA seq standardized for log
                      type_key = 'celltype',        # Keywords representing cell types (in protein dataset)
                      data_dir="/scratch/jfisher2/temp/scMMT_preprocessed.pkl",  # Save path for processed data
                      data_load=False,                # Do you want to import existing processed data
                      dataset_batch = True,           # Is there a batch effect in the training set and testing machine
                      log_weight=3,                   # Log weights for different cell types
                      val_split = None,               # Do you need to divide the validation set according to the distribution of the test set
                      min_cells = 0,                  # Minimum cell count filtering
                      min_genes = 0,                  # Minimum number of genes filtering
                      n_svd = 300,                    # Dimension obtained using Tsvd dimensionality reduction
                      n_fa=180,                       # Dimension obtained by using FA dimensionality reduction
                      n_hvg=550,                      # Number of high variants obtained through screening
                     )



scMMT.train(n_epochs = 100, ES_max = 12, decay_max = 6, decay_step = 0.1, lr = 10**(-3), label_smoothing=0.4, 
            h_size=600, drop_rate=0.15, n_layer=4,
            weights_dir = "model_weight", load = False)


predicted_test = scMMT.predict()


predicted_test.write_h5ad(
    args.basedir + "/output/scMMT/scMMT_prediction.h5ad"
)

if dobenchmark=='true':
    #get benchmarking prefix
    #make sure looking at filename only, no extended file path
    basename=os.path.basename(rna_train_data_file)
    prefix= basename.replace("_training_data_rna_norm.csv", "") 
    np.savetxt(args.basedir + "/output/scMMT/" +prefix + "scMMT_prediction.csv", predicted_test.X, delimiter=",")
    np.savetxt(prefix + "scMMT_prediction.csv", predicted_test.X, delimiter=",")
else:
    #save for later
    np.savetxt(args.basedir + "/output/scMMT/scMMT_prediction.csv", predicted_test.X, delimiter=",")

    #pass to pipeline
    np.savetxt("scMMT_prediction.csv", predicted_test.X, delimiter=",")
