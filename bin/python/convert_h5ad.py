import numpy as np
import scanpy as sc
import scipy.sparse as sp
import argparse
import anndata
import pandas as pd
from anndata import AnnData, read_h5ad
from scanpy import read
import os




#setup arg parse
parser = argparse.ArgumentParser( prog = 'Script to train scipenn on scRNAseq data')

parser.add_argument('-d', '--basedir', required=True, help="pipeline base directory")
parser.add_argument('-f','--files', nargs='+', help='<Required> Set flag', required=True)


args = parser.parse_args()

print(args)

crossmodalnet_dir=args.basedir +  '/output/CrossModalNet/cite_touse2'
if not os.path.exists(crossmodalnet_dir):
    os.makedirs(crossmodalnet_dir)


input_files = args.files

#need to clean inputs of wrong character
characters_to_remove = ['[', ']', ',']
translation_table = str.maketrans('', '', ''.join(characters_to_remove))

input_files = args.files

metadata_file =  [file for file in input_files if 'metadata' in file][0]
metadata_file = metadata_file.translate(translation_table)

rna_train_data_file = [file for file in input_files if 'training_data_rna_raw' in file][0]
rna_train_data_file = rna_train_data_file.translate(translation_table)

prot_train_data_file = [file for file in input_files if 'training_data_prot_raw' in file][0]
prot_train_data_file = prot_train_data_file.translate(translation_table)

rna_test_data_file =  [file for file in input_files if 'testing_data_rna_raw' in file][0]
rna_test_data_file = rna_test_data_file.translate(translation_table)





adata_rna_train = anndata.read_csv(rna_train_data_file, delimiter=',', first_column_names = True).transpose()
adata_prot_train = anndata.read_csv(prot_train_data_file, delimiter=',', first_column_names = True).transpose()
adata_rna_test = anndata.read_csv(rna_test_data_file, delimiter=',', first_column_names = True).transpose()

# # # read in as pandas
meta = pd.read_csv(metadata_file, index_col=0)


#add metadata
ind=meta.index.isin(adata_rna_train.obs.index)
adata_rna_train.obs['donor']=meta.loc[ind,'donor']
adata_rna_train.obs['celltype']=meta.loc[ind,'celltype']

ind=meta.index.isin(adata_prot_train.obs.index)
adata_prot_train.obs['donor']=meta.loc[ind,'donor']
adata_prot_train.obs['celltype']=meta.loc[ind,'celltype']

ind=meta.index.isin(adata_rna_test.obs.index)
adata_rna_test.obs['donor']=meta.loc[ind,'donor']
adata_rna_test.obs['celltype']=meta.loc[ind,'celltype']

#fix metadata, convert to categorical dtype
adata_rna_train.obs['donor'] = adata_rna_train.obs['donor'].astype('category')
adata_prot_train.obs['donor'] = adata_prot_train.obs['donor'].astype('category')
adata_rna_test.obs['donor'] = adata_rna_test.obs['donor'].astype('category')
adata_rna_train.obs['celltype'] = adata_rna_train.obs['celltype'].astype('category')
adata_prot_train.obs['celltype'] = adata_prot_train.obs['celltype'].astype('category')
adata_rna_test.obs['celltype'] = adata_rna_test.obs['celltype'].astype('category')

#add timepoint category for CrossModalNet to check
adata_rna_train.obs['day']=0
adata_prot_train.obs['day'] = 0
adata_rna_test.obs['day'] = 0


adata_rna_train.var['features'] = pd.DataFrame(adata_rna_train.var.index, index=adata_rna_train.var.index)
adata_rna_train.var['features'] = adata_rna_train.var['features'].astype('category')
adata_prot_train.var['features'] = pd.DataFrame(adata_prot_train.var.index, index=adata_prot_train.var.index)
adata_prot_train.var['features'] = adata_prot_train.var['features'].astype('category')
adata_rna_test.var['features'] = pd.DataFrame(adata_rna_test.var.index, index=adata_rna_test.var.index)
adata_rna_test.var['features'] = adata_rna_test.var['features'].astype('category')




#save as h5ad and output
adata_rna_train.write_h5ad(
    args.basedir +  '/output/training_files/training_rna_data.h5ad'
)

adata_prot_train.write_h5ad(
    args.basedir +  '/output/training_files/training_prot_data.h5ad'
)

adata_rna_test.write_h5ad(
    args.basedir +  '/output/training_files/testing_rna_data.h5ad'
)





adata_rna_train.write_h5ad(
    'training_rna_data.h5ad'
)

adata_prot_train.write_h5ad(
    'training_prot_data.h5ad'
)

adata_rna_test.write_h5ad(
    'testing_rna_data.h5ad'
)
