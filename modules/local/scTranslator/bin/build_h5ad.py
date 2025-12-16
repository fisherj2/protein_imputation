

import warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import pickle

import anndata as ad
import argparse
import os


parser = argparse.ArgumentParser( prog = 'Script to  build h5ad ready for scTranslator')
parser.add_argument('-d', '--basedir', required=True, help="pipeline base directory")
parser.add_argument('-b','--bench',  help='<Required> Set flag for benchmarking', required=True)
parser.add_argument('-f','--files', nargs='+', help='<Required> Set flag', required=True)
parser.add_argument('-m','--moduledir', nargs='+', help='directory of the scTranslator module', required=True)     
args = parser.parse_args()
warnings.filterwarnings('ignore')
    
dobenchmark = args.bench
input_files = args.files
basedir = args.basedir
moduleDir = args.moduledir[0]

#if length is one, probably need to split
if len(input_files)==1 or isinstance(input_files, (str)):
    print('adjusting input format')
    if isinstance(input_files, (list)):
        input_files=input_files[0]
    
    input_files=input_files.split('.csv')
    input_files=[s + '.csv' for s in input_files]

print(input_files)

# #check output dir exists
# if not os.path.exists(args.launchdir + "/output/scTranslator"):
#     os.makedirs(args.launchdir + "/output/scTranslator")

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

prot_test_data_file =  [file for file in input_files if 'testing_data_prot_raw' in file][0]
prot_test_data_file = prot_test_data_file.translate(translation_table)

test_rna_var_meta_file =  [file for file in input_files if 'testing_data_RNA_var_meta' in file][0]
test_rna_var_meta_file = test_rna_var_meta_file.translate(translation_table)

test_prot_var_meta_file =  [file for file in input_files if 'testing_data_protein_var_meta' in file][0]
test_prot_var_meta_file = test_prot_var_meta_file.translate(translation_table)

train_rna_var_meta_file =  [file for file in input_files if 'training_data_RNA_var_meta' in file][0]
train_rna_var_meta_file = train_rna_var_meta_file.translate(translation_table)

train_prot_var_meta_file =  [file for file in input_files if 'training_data_protein_var_meta' in file][0]
train_prot_var_meta_file = train_prot_var_meta_file.translate(translation_table)

#remove whitespace
metadata_file_train=metadata_file_train.strip()
metadata_file_test=metadata_file_test.strip()
rna_train_data_file=rna_train_data_file.strip()
prot_train_data_file=prot_train_data_file.strip()
rna_test_data_file=rna_test_data_file.strip()
prot_test_data_file=prot_test_data_file.strip()

test_prot_var_meta_file=test_prot_var_meta_file.strip()
test_rna_var_meta_file=test_rna_var_meta_file.strip()
train_prot_var_meta_file=train_prot_var_meta_file.strip()
train_rna_var_meta_file=train_rna_var_meta_file.strip()

#normalise ready for scTranslator

# metadata_file =  '/shared-workspace/jfisher2/analyses/protein_prediction/testing/output/training_files/0.5_1_metadata_train.csv'
# rna_train_data_file = '/shared-workspace/jfisher2/analyses/protein_prediction_publication/testing/output/training_files/0.5_1_training_data_rna_raw.csv'
# prot_train_data_file = '/shared-workspace/jfisher2/analyses/protein_prediction_publication/testing/output/training_files/0.5_1_training_data_prot_raw.csv'
# rna_test_data_file =  '/shared-workspace/jfisher2/analyses/protein_prediction_publication/testing/output/testing_files/0.5_1_testing_data_rna_raw.csv'
# prot_test_data_file = '/shared-workspace/jfisher2/analyses/protein_prediction_publication/testing/output/testing_files/0.5_1_testing_data_prot_raw.csv'

train_rna = ad.read_csv(rna_train_data_file, delimiter=',', first_column_names = True).transpose()
train_protein = ad.read_csv(prot_train_data_file, delimiter=',', first_column_names = True).transpose()
test_rna = ad.read_csv(rna_test_data_file, delimiter=',', first_column_names = True).transpose()
test_protein = ad.read_csv(prot_test_data_file, delimiter=',', first_column_names = True).transpose()

test_prot_var_meta = pd.read_csv(test_prot_var_meta_file)
test_rna_var_meta = pd.read_csv(test_rna_var_meta_file)
train_prot_var_meta = pd.read_csv(train_prot_var_meta_file)
train_rna_var_meta = pd.read_csv(train_rna_var_meta_file)

#check if 'scTranslator_ID' exists as a column
if 'scTranslator_ID' not in test_prot_var_meta.columns:
    print("Error: Column 'scTranslator_ID' not found! Make sure it exists in both your RNA and protein feature metadata, see scTranslator github for ID contents.")
    sys.exit(1)
    
if 'scTranslator_ID' not in test_rna_var_meta.columns:
    print("Error: Column 'scTranslator_ID' not found! Make sure it exists in both your RNA and protein feature metadata, see scTranslator github for ID contents.")
    sys.exit(1)
    
test_protein.var['my_Id'] = test_prot_var_meta.set_index(test_prot_var_meta.columns[0]).loc[test_protein.var.index, 'scTranslator_ID']
#test_protein.var['my_Id'] = test_prot_var_meta['scTranslator_ID'].values
#test_protein.var['my_Id'] = test_protein.var['my_Id'].astype('str')
test_rna.var['my_Id'] = test_rna_var_meta.set_index(test_rna_var_meta.columns[0]).loc[test_rna.var.index, 'scTranslator_ID']
#test_rna.var['my_Id'] = test_rna_var_meta['scTranslator_ID'].values
#test_rna.var['my_Id'] = test_rna.var['my_Id'].astype('str')
train_protein.var['my_Id'] = train_prot_var_meta.set_index(train_prot_var_meta.columns[0]).loc[train_protein.var.index, 'scTranslator_ID']
#train_protein.var['my_Id'] = train_prot_var_meta['scTranslator_ID'].values
#train_protein.var['my_Id'] = train_protein.var['my_Id'].astype('str')
train_rna.var['my_Id'] = train_rna_var_meta.set_index(train_rna_var_meta.columns[0]).loc[train_rna.var.index, 'scTranslator_ID']
#train_rna.var['my_Id'] = train_rna_var_meta['scTranslator_ID'].values
#train_rna.var['my_Id'] = train_rna.var['my_Id'].astype('str')

test_protein.var['orig_feature'] = test_protein.var.index
test_rna.var['orig_feature'] = test_rna.var.index
train_protein.var['orig_feature'] = train_protein.var.index
train_rna.var['orig_feature'] = train_rna.var.index


test_protein = test_protein[:, test_protein.var['my_Id'].notna()]
test_rna = test_rna[:, test_rna.var['my_Id'].notna()]
train_protein = train_protein[:, train_protein.var['my_Id'].notna()]
train_rna = train_rna[:, train_rna.var['my_Id'].notna()]


# test_protein.var.index = test_protein.var['my_Id']
# test_rna.var.index = test_rna.var['my_Id']
# train_protein.var.index = train_protein.var['my_Id']
# train_rna.var.index = train_rna.var['my_Id']
# 


# #set all values of the test object to zero, this will be a dummy object for prediction
# test_protein.X = np.zeros_like(test_protein.X)
# 
# #some proteins map to the same gene symbol. In this case remove the duplicates. 
# test_protein = test_protein[:, ~test_protein.var_names.duplicated(keep='first')]

# # Create protein query the way scTranslator expects
# obs = pd.DataFrame(test_rna.obs.values.tolist(), index=test_rna.obs.index)
# 
# 
# #remove duplicate features
# is_first_occurrence = ~test_protein.var_names.duplicated(keep='first')
# 
# # Subset the AnnData object
# test_protein = test_protein[:, is_first_occurrence].copy()
# 
# is_first_occurrence = ~train_protein.var_names.duplicated(keep='first')
# 
# # Subset the AnnData object
# train_protein = train_protein[:, is_first_occurrence].copy()

#test_protein.X = np.zeros((test_protein.n_obs, test_protein.var.shape[0]))


#in the case of the training protein data we don't need to remove the duplicates

#train_protein = train_protein[:, ~train_protein.var.index.duplicated(keep=False)]

if dobenchmark=='true':
    #get benchmarking prefix
    #make sure looking at filename only, no extended file path
    basename=os.path.basename(rna_train_data_file)
    prefix= basename.replace("_training_data_rna_raw.csv", "") 
    #write to h5ad
    train_rna.write(prefix +'_rna_train_raw.h5ad')
    train_protein.write(prefix +'_protein_train_raw.h5ad')
    test_rna.write(prefix+'_rna_test_raw.h5ad')
    test_protein.write(prefix+'_protein_test_raw.h5ad')
    #save for review
    # train_rna.write(args.launchdir + "/output/training_files/" + prefix + "rna_train_raw.h5ad")
    # train_protein.write(args.launchdir + "/output/training_files/" + prefix + "protein_train_raw.h5ad")
    # test_rna.write(args.launchdir + "/output/training_files/" + prefix + "rna_test_raw.h5ad")
else:
   #write to h5ad
    train_rna.write('rna_train_raw.h5ad')
    train_protein.write('protein_train_raw.h5ad')
    test_rna.write('rna_test_raw.h5ad')
    test_protein.write('protein_test_raw.h5ad')
    #save for review
    # train_rna.write(args.launchdir + "/output/training_files/rna_train_raw.h5ad")
    # train_protein.write(args.launchdir + "/output/training_files/protein_train_raw.h5ad")
    # test_rna.write(args.launchdir + "/output/training_files/rna_test_raw.h5ad")
    # 
    # 

