import numpy as np
import torch
import scanpy as sc
import scvi
from scvi.model import TOTALVI

import scipy.sparse as sp
import argparse
import anndata as ad
import pandas as pd
from anndata import AnnData, read_h5ad
from scanpy import read
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



#CHECK OUTPUT DIR EXISTS
model_directory=args.launchdir + '/output/totalVI'
if not os.path.exists(model_directory):
    os.makedirs(model_directory)



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

# #convert to sparse matrix
# adata_rna_train.X = sp.csr_matrix(adata_rna_train.X)
# adata_prot_train.X = sp.csr_matrix(adata_prot_train.X)
# adata_rna_test.X = sp.csr_matrix(adata_rna_test.X)


#create one training object with both rna and protein data
adata_train=adata_rna_train
#create as a pandas dataframe
adata_train.obsm["protein_expression"]=pd.DataFrame(adata_prot_train.X, columns = adata_prot_train.var['features'].index, index=adata_prot_train.obs_names)

model_directory = args.launchdir + "/output/totalVI"




ref = adata_train
query = adata_rna_test

query.obs["batch"] = "query_data"
ref.obs["batch"] = "ref_data"

# put matrix of zeros for protein expression (considered missing)
pro_exp = ref.obsm["protein_expression"]
data = np.zeros((query.n_obs, pro_exp.shape[1]))
query.obsm["protein_expression"] = pd.DataFrame(
    index=query.obs_names, data=data, columns=pro_exp.columns
)



sc.pp.highly_variable_genes(
    ref,
    n_top_genes=4000,
    flavor="seurat_v3",
    batch_key="batch",
    subset=True,
)

query = query[:, ref.var_names].copy()


scvi.model.TOTALVI.setup_anndata(
    ref, batch_key="batch", protein_expression_obsm_key="protein_expression"
)

totalvi_ref = scvi.model.TOTALVI(ref, use_layer_norm="both", use_batch_norm="none")

totalvi_ref.train()


TOTALVI_LATENT_KEY = "X_totalVI"

ref.obsm[TOTALVI_LATENT_KEY] = totalvi_ref.get_latent_representation()

totalvi_ref.save(model_directory, overwrite=True)


scvi.model.TOTALVI.prepare_query_anndata(query, model_directory)
totalvi_query = scvi.model.TOTALVI.load_query_data(
    query,
    model_directory,
)

totalvi_query.train(200, plan_kwargs={"weight_decay": 0.0})

query.obsm[TOTALVI_LATENT_KEY] = totalvi_query.get_latent_representation()

_, imputed_proteins = totalvi_query.get_normalized_expression(
    query,
    n_samples=10,
    return_mean=True,
    transform_batch=["ref_data", "query_data"],
)



if dobenchmark=='true':
    #get benchmarking prefix
    #make sure looking at filename only, no extended file path
    basename=os.path.basename(rna_train_data_file)
    prefix= basename.replace("_training_data_rna_norm.csv", "") 
    np.savetxt(args.launchdir + "/output/totalVI/" +prefix + "totalVI_prediction.csv", imputed_proteins, delimiter=",")
    np.savetxt(prefix + "totalVI_prediction.csv", imputed_proteins, delimiter=",")
else:
    #save for later
    np.savetxt(args.launchdir + "/output/totalVI/totalVI_prediction.csv", imputed_proteins, delimiter=",")

    #pass to pipeline
    np.savetxt("totalVI_prediction.csv", imputed_proteins, delimiter=",")
