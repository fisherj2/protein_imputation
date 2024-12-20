
print('loading libraries')
import argparse
import pandas as pd
import numpy as np
import os
import scanpy as sc
import numpy as np
import scipy.sparse as sp
import anndata as ad
from evaluate import evaluate
from prediction import ADTPredictor
print('libraries imported')

print('parsing inputs')
#setup arg parse
parser = argparse.ArgumentParser( prog = 'Script to train sclinear on scRNAseq data')
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
if not os.path.exists(args.basedir + "/output/scLinear"):
    os.makedirs(args.basedir + "/output/scLinear")

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

#convert to sparse matrix
adata_rna_train.X = sp.csc_matrix(adata_rna_train.X)
adata_prot_train.X = sp.csc_matrix(adata_prot_train.X)
adata_rna_test.X = sp.csc_matrix(adata_rna_test.X)

print('running prediction')
pipe = ADTPredictor(do_log1p=True)
pipe.fit(gex_train=adata_rna_train.X.toarray(), adt_train=adata_prot_train.X.toarray(), gex_test=adata_rna_test.X.toarray(), gex_names=adata_rna_train.var_names.to_numpy(), adt_names=adata_prot_train.var_names.to_numpy())
adt_pred, adt_names = pipe.predict(adata_rna_test.X.toarray())

#also compute feature importance. Turned off here due to large memory consumption
#impvals = pipe.feature_importance( gex_test= adata_rna_test.X)


if dobenchmark=='true':
    #get benchmarking prefix
    #make sure looking at filename only, no extended file path
    basename=os.path.basename(rna_train_data_file)
    prefix= basename.replace("training_data_rna_raw.csv", "") 
    np.savetxt(args.basedir + "/output/scLinear/" +prefix + "sclinear_prediction.csv", adt_pred, delimiter=",")
    np.savetxt(prefix + "sclinear_prediction.csv", adt_pred, delimiter=",")
    pipe.save(path =args.basedir + "/output/scLinear/" +prefix + "sclinear_model" )
    #np.savetxt(args.basedir + "/output/scLinear/" +prefix + "sclinear_model_feature_importance.csv" ,impvals,delimiter=",")
else:
    #save for later
    np.savetxt(args.basedir + "/output/scLinear/sclinear_prediction.csv", adt_pred, delimiter=",")
    #np.savetxt(args.basedir + "/output/scLinear/sclinear_model_feature_importance.csv" ,impvals,delimiter=",")

    #pass to pipeline
    np.savetxt("sclinear_prediction.csv", adt_pred, delimiter=",")
    pipe.save(path =args.basedir + "/output/scLinear/sclinear_model" )
