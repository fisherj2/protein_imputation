
import os
import time
import datetime
import argparse
import warnings

import anndata as ad
import torch
import torch.optim as optim
from torch.optim.lr_scheduler import StepLR


import torch.optim as optim
import scanpy as sc
import numpy as np
import pandas as pd
from sklearn.model_selection import ShuffleSplit

import sys 

#set seed
seed = 123
torch.manual_seed(seed)
torch.cuda.manual_seed(seed) 
torch.cuda.manual_seed_all(seed)  # if you are using multi-GPU
np.random.seed(seed)

parser = argparse.ArgumentParser(description='PyTorch Example')
parser.add_argument('--repeat', type=int, default=1,
                    help='for repeating experiments to change seed (default: 1)')
parser.add_argument('--test_batch_size', type=int, default=4,
                    help='input batch size for testing (default: 4)')
parser.add_argument('--no-cuda', action='store_true', default=False,
                    help='disables CUDA training')
parser.add_argument('--seed', type=int, default=1105,
                    help='random seed (default: 1105)')
parser.add_argument('--enc_max_seq_len', type=int, default=20000,
                    help='sequence length of encoder')
parser.add_argument('--dec_max_seq_len', type=int, default=1000,
                    help='sequence length of decoder')
parser.add_argument('--fix_set', action='store_false',
                    help='fix (aligned) or disordering (un-aligned) dataset')
parser.add_argument('--pretrain_checkpoint', default='checkpoint/stage2_single-cell_scTranslator.pt',
                    help='path for loading the pretrain checkpoint')
parser.add_argument('--RNA_path', default='dataset/test/dataset1/GSM5008737_RNA_finetune_withcelltype.h5ad',
                    help='path for loading the rna')
parser.add_argument('--Pro_path', default='dataset/test/dataset1/GSM5008738_protein_finetune_withcelltype.h5ad',
                    help='path for loading the protein')
parser.add_argument('-d', '--basedir', required=True, help="pipeline base directory")
parser.add_argument('-b','--bench',  help=' Set flag for benchmarking', required=True)
parser.add_argument('-f','--files', nargs='+', help='input files', required=True)                   
parser.add_argument('-m','--moduledir', nargs='+', help='directory of the scTranslator module', required=True)     
args = parser.parse_args()
warnings.filterwarnings('ignore')

dobenchmark = args.bench
input_files = args.files
moduleDir = args.moduledir[0]
# 
# #check output dir exists
# if not os.path.exists(args.launchdir + "/output/scTranslator"):
#     os.makedirs(args.launchdir + "/output/scTranslator")
#     

# Set the correct path
model_path = os.path.join(moduleDir+ '/bin/code/model')
print(model_path)

# Add to Python path
sys.path.insert(0, model_path)

from performer_enc_dec import *
from utils import *


###########################
#--- Prepare The Model ---#
###########################
device = 'cuda' if torch.cuda.is_available() else 'cpu'
print('device',device)
model = torch.load(args.pretrain_checkpoint, map_location=torch.device(device))
# model = model.to(device)





#if length is one, probably need to split
if len(input_files)==1 or isinstance(input_files, (str)):
    print('adjusting input format')
    if isinstance(input_files, (list)):
        input_files=input_files[0]
    
    input_files=input_files.split('.h5ad')
    input_files=[s + '.h5ad' for s in input_files]

input_files = [s.strip() for s in input_files]
print(input_files)

##########################
#--- Prepare The Data ---#
##########################
     
#---  Load Single Cell Data  ---#
# scRNA_adata = sc.read_h5ad(args.RNA_path)[:100]
# scP_adata = sc.read_h5ad(args.Pro_path)[:100]
# print('Total number of origin RNA genes: ', scRNA_adata.n_vars)
# print('Total number of origin proteins: ', scP_adata.n_vars)
# print('Total number of origin cells: ', scRNA_adata.n_obs)
# print('# of NAN in X', np.isnan(scRNA_adata.X).sum())
# print('# of NAN in X', np.isnan(scP_adata.X).sum())

#need to clean inputs of wrong character
characters_to_remove = ['[', ']', ',']
translation_table = str.maketrans('', '', ''.join(characters_to_remove))


# rna_train_data_file = [file for file in input_files if 'rna_train_raw' in file][0]
# rna_train_data_file = rna_train_data_file.translate(translation_table)
# 
# protein_train_data_file = [file for file in input_files if 'protein_train_raw' in file][0]
# protein_train_data_file = protein_train_data_file.translate(translation_table)

rna_test_data_file =  [file for file in input_files if 'rna_test_raw' in file][0]
rna_test_data_file = rna_test_data_file.translate(translation_table).strip() 

protein_test_data_file =  [file for file in input_files if 'protein_test_raw' in file][0]
protein_test_data_file = protein_test_data_file.translate(translation_table).strip() 



test_rna = sc.read_h5ad(rna_test_data_file)

test_protein = sc.read_h5ad(protein_test_data_file)

def find_optimal_batch_size(num_samples, preferred_batch_size, min_remainder=2):
    """Find batch size that avoids single-sample batches"""
    
    # First try the preferred batch size
    remainder = num_samples % preferred_batch_size
    if remainder == 0:
        return preferred_batch_size
    
    # Find a good divisor
    for batch_size in range(10, 1, -1):
        remainder = num_samples % batch_size
        
        if remainder == 0 :
            return batch_size
    
    # Fallback
    print('no better batch size found, returning default')
    return preferred_batch_size

# Get the actual number of cells/samples
num_test_samples = test_protein.n_obs  # Number of cells, not features
optimal_batch_size = find_optimal_batch_size(num_test_samples, args.test_batch_size,min_remainder=args.test_batch_size)

if num_test_samples % optimal_batch_size == 1:
    optimal_batch_size += 1
    
print("using batch size " +  str(optimal_batch_size))

# #---  Construct Dataloader ---#

if args.fix_set == True:
    my_testset = fix_SCDataset(test_rna, test_protein, args.enc_max_seq_len, args.dec_max_seq_len)
else:
    my_testset = SCDataset(test_rna, test_protein, args.enc_max_seq_len, args.dec_max_seq_len)


my_testset = SCDataset(test_rna, test_protein, args.enc_max_seq_len, args.dec_max_seq_len)

test_loader = torch.utils.data.DataLoader(my_testset, batch_size=optimal_batch_size, drop_last=False)
print("load data ended")

##################
#---  Testing ---#
##################
start_time = time.time()


test_loss, test_ccc, test_pcor, y_hat, y, indices = test(model, device, test_loader)

test_cells = test_protein.obs.index
test_feat = test_protein.var['orig_feature']

y_pred =  pd.DataFrame(y_hat,columns=test_feat, index=test_cells)
y_truth = pd.DataFrame(y,columns=test_feat, index=test_cells)
indices = pd.DataFrame(indices)

#set rownames and colnames
test_pcor = pd.DataFrame({
    'feature': test_feat,
    'pearson_correlation': test_pcor
}).set_index('feature')
# test_loss = pd.DataFrame({
#     'feature': test_feat,
#     'pearson_correlation': test_loss
# }).set_index('feature')
# 
# test_ccc = pd.DataFrame({
#     'feature': test_feat,
#     'pearson_correlation': test_ccc
# }).set_index('feature')
# #set rownames and colnames


##############################
#---  Prepare for Storage ---#
##############################


if args.RNA_path == 'dataset/test/dataset1/GSM5008737_RNA_finetune_withcelltype.h5ad':
    dataset_flag = '/seuratv4_16W-without_fine-tune'
else:
    dataset_flag = '/new_data-without_fine-tune'
file_path = 'result/test'+dataset_flag
if not os.path.exists(file_path):
    os.makedirs(file_path)

dict = vars(args)
filename = open(file_path+'/args'+str(args.repeat)+'.txt','w')
for k,v in dict.items():
    filename.write(k+':'+str(v))
    filename.write('\n')
filename.close()

#---  Save the Final Results ---#
# log_path = file_path+'/performance_log.csv'
# log_all = pd.DataFrame(columns=['test_loss', 'test_ccc'])
# log_all.loc[args.repeat] = np.array([test_loss, test_ccc])
# log_all.to_csv(log_path)
# y_pred.to_csv(file_path+'/y_pred.csv')
# y_truth.to_csv(file_path+'/y_truth.csv')



if dobenchmark=='true':
    #get benchmarking prefix
    #make sure looking at filename only, no extended file path
    basename=os.path.basename(rna_test_data_file)
    prefix= basename.replace("_testing_data_rna_raw.csv", "") 
    
    #log_path = args.launchdir + '/output/scTranslator/' + prefix + '_nofinetune_performance_log.csv'
    log_path = prefix + '_nofinetune_performance_log.csv'
    log_all = pd.DataFrame(columns=['test_loss', 'test_ccc'])
    log_all.loc[args.repeat] = np.array([test_loss, test_ccc])
    log_all.to_csv(log_path)

    test_pcor.to_csv(prefix + '_nofinetune_pearson_corr.csv')
    # test_ccc.to_csv(prefix + '_nofinetune_ccc.csv')
    # test_loss.to_csv(prefix + '_nofinetune_loss.csv')

    y_pred.to_csv( prefix + "scTranslator_nofinetune_prediction.csv")
    y_truth.to_csv(prefix + "scTranslator_nofinetune_truth.csv")
    indices.to_csv(prefix + "scTranslator_cell_indices.csv")
else:
  
    #log_path = args.launchdir + '/output/scTranslator/' + '_nofinteune_performance_log.csv'
    log_path =  'nofinteune_performance_log.csv'
    log_all = pd.DataFrame(columns=['test_loss', 'test_ccc'])
    log_all.loc[args.repeat] = np.array([test_loss, test_ccc])
    log_all.to_csv(log_path)
    test_pcor.to_csv('nofinetune_pearson_corr.csv')
    # test_ccc.to_csv( 'nofinetune_ccc.csv')
    # test_loss.to_csv('nofinetune_loss.csv')

    
    #pass downstream to pipe
    y_pred.to_csv( "scTranslator_nofinetune_prediction.csv")
    y_truth.to_csv("scTranslator_nofinetune_truth.csv")
    indices.to_csv("scTranslator_cell_indices.csv")
    
