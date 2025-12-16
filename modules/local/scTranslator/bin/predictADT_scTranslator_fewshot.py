
import os
import time
import datetime
import argparse
import warnings
import hashlib

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
parser.add_argument('--batch_size', type=int, default=1, metavar='N',
                    help='input batch size for each GPU training (default: 1)')
parser.add_argument('--test_batch_size', type=int, default=4,
                    help='input batch size for testing (default: 4)')
parser.add_argument('--epochs', type=int, default=100, metavar='N',
                    help='number of epochs to train (default: 100)')
parser.add_argument('--lr', type=float, default=2*1e-4, metavar='LR',
                    help='learning rate (default: 1.0)')
parser.add_argument('--gamma', type=float, default=1, metavar='M',
                    help='Learning rate step gamma (default: 1 (not used))')
parser.add_argument('--gamma_step', type=float, default=2000,
                    help='Learning rate step (default: 2000 (not used))')
parser.add_argument('--no-cuda', action='store_true', default=False,
                    help='disables CUDA training')
parser.add_argument('--seed', type=int, default=1105,
                    help='random seed (default: 1105)')
parser.add_argument('--repeat', type=int, default=1,
                    help='for repeating experiments to change seed (default: 1)')
parser.add_argument('--local_rank', default=0, type=int,
                    help='node rank for distributed training')
parser.add_argument('--frac_finetune_test', type=float, default=0.1,
                    help='test set ratio')
parser.add_argument('--dim', type=int, default=128,
                    help='latend dimension of each token')
parser.add_argument('--enc_max_seq_len', type=int, default=20000,
                    help='sequence length of encoder')
parser.add_argument('--dec_max_seq_len', type=int, default=1000,
                    help='sequence length of decoder')
parser.add_argument('--translator_depth', type=int, default=2,
                    help='translator depth')
parser.add_argument('--initial_dropout', type=float, default=0.1,
                    help='sequence length of decoder')
parser.add_argument('--enc_depth', type=int, default=2,
                    help='sequence length of decoder')
parser.add_argument('--enc_heads', type=int, default=8,
                    help='sequence length of decoder')
parser.add_argument('--dec_depth', type=int, default=2,
                    help='sequence length of decoder')
parser.add_argument('--dec_heads', type=int, default=8,
                    help='sequence length of decoder')
parser.add_argument('--fix_set', action='store_false',
                    help='fix (aligned) or disordering (un-aligned) dataset')
parser.add_argument('--pretrain_checkpoint', default='checkpoint/stage2_single-cell_scTranslator.pt',
                    help='path for loading the pretrain checkpoint')
parser.add_argument('--resume', default=False, help='resume training from breakpoint')
parser.add_argument('--path_checkpoint', default='checkpoint/stage2_single-cell_scTranslator.pt',
                    help='path for loading the resume checkpoint (need specify)')
parser.add_argument('--RNA_path', default='dataset/test/dataset1/GSM5008737_RNA_finetune_withcelltype.h5ad',
                    help='path for loading the rna')
parser.add_argument('--Pro_path', default='dataset/test/dataset1/GSM5008738_protein_finetune_withcelltype.h5ad',
                    help='path for loading the protein')
parser.add_argument('-d', '--basedir', required=True, help="pipeline base directory")
parser.add_argument('-b','--bench',  help=' Set flag for benchmarking', required=True)
parser.add_argument('-f','--files', nargs='+', help='input files', required=True)                   
parser.add_argument('-m','--moduledir', nargs='+', help='directory of the scTranslator module', required=True)                   
parser.add_argument('--ncells', default=100,
                    help='number of cells to use for finetuning')
args = parser.parse_args()
 
dobenchmark = args.bench
input_files = args.files

#if length is one, probably need to split
if len(input_files)==1 or isinstance(input_files, (str)):
    print('adjusting input format')
    if isinstance(input_files, (list)):
        input_files=input_files[0]
    
    input_files=input_files.split('.h5ad')
    input_files=[s + '.h5ad' for s in input_files]

input_files = [s.strip() for s in input_files]
print(input_files)

moduleDir = args.moduledir[0]

# #check output dir exists
# if not os.path.exists(args.launchdir + "/output/scTranslator"):
#     os.makedirs(args.launchdir + "/output/scTranslator")
    
# Set the correct path
model_path = os.path.expanduser(moduleDir+ '/bin/code/model')

# Add to Python path
sys.path.insert(0, model_path)

from performer_enc_dec import *
from utils import *

#########################
#--- Prepare for DDP ---#
#########################
use_cuda = not args.no_cuda and torch.cuda.is_available()
print("use_cuda: %s" % use_cuda)
ngpus_per_node = torch.cuda.device_count()
print("ngpus_per_node: %s" % ngpus_per_node)
is_distributed = ngpus_per_node > 1
print('seed', args.seed)
setup_seed(args.seed)
print(torch.__version__)
# Initializes the distributed environment to help process communication
torch.distributed.init_process_group(backend='nccl', timeout=datetime.timedelta(seconds=15400))
# Each process sets the GPU it should use based on its local rank
print("local_rank: %s" % args.local_rank)
device = torch.device("cuda", args.local_rank)
print(device)
torch.cuda.set_device(args.local_rank)
rank = int(os.environ['RANK'])
print('rank', rank)

###########################
#--- Prepare The Model ---#
###########################
model = scPerformerEncDec(
    dim=args.dim,
    translator_depth=args.translator_depth,
    initial_dropout=args.initial_dropout,
    enc_depth=args.enc_depth,
    enc_heads=args.enc_heads,
    enc_max_seq_len=args.enc_max_seq_len,
    dec_depth=args.dec_depth,
    dec_heads=args.dec_heads,
    dec_max_seq_len=args.dec_max_seq_len
    )
model = torch.load(args.pretrain_checkpoint)
# Resume training from breakpoints
if args.resume == True:
    checkpoint = torch.load(args.path_checkpoint)
    model = checkpoint['net']
    model = model.to(device)
    if is_distributed:
        print("start init process group")
        # device_ids will include all GPU devices by default
        model = torch.nn.parallel.DistributedDataParallel(model, device_ids=[args.local_rank], find_unused_parameters=True)
        print("end init process group")
    #---  Prepare Optimizer ---#
    optimizer = optim.Adam(model.parameters(), lr=args.lr, amsgrad=True) 
    optimizer.load_state_dict(checkpoint['optimizer'])
    #---  Prepare Scheduler ---#
    scheduler = StepLR(optimizer, step_size=args.gamma_step, gamma=args.gamma)
    scheduler.load_state_dict(checkpoint['scheduler'])
    start_epoch = checkpoint['epoch']
else:
    start_epoch = 0
    model = model.to(device)
    if is_distributed:
        print("start init process group")
        # device_ids will include all GPU devices by default
        model = torch.nn.parallel.DistributedDataParallel(model, device_ids=[args.local_rank], find_unused_parameters=True)
        print("end init process group")
    #---  Prepare Optimizer ---#
    optimizer = optim.Adam(model.parameters(), lr=args.lr, amsgrad=True)
    #---  Prepare Scheduler ---#
    scheduler = StepLR(optimizer, step_size=args.gamma_step, gamma=args.gamma)
    
    
        
        
        
##########################
#--- Prepare The Data ---#
##########################
     
# #---  Load Single Cell Data  ---#
# scRNA_adata = sc.read_h5ad(args.RNA_path)[:100]
# scP_adata = sc.read_h5ad(args.Pro_path)[:100]
# print('Total number of origin RNA genes: ', scRNA_adata.n_vars)
# print('Total number of origin proteins: ', scP_adata.n_vars)
# print('Total number of origin cells: ', scRNA_adata.n_obs)
# print('# of NAN in X', np.isnan(scRNA_adata.X).sum())
# print('# of NAN in X', np.isnan(scP_adata.X).sum())
# 
# #---  Seperate Training and Testing set ---#
# setup_seed(args.seed+args.repeat)
# train_index, test_index = next(ShuffleSplit(n_splits=1,test_size=args.frac_finetune_test).split(scRNA_adata.obs.index))
# # --- RNA ---#
# train_rna = scRNA_adata[train_index]
# test_rna = scRNA_adata[test_index]
# # --- Protein ---#
# train_protein = scP_adata[train_index]
# test_protein = scP_adata[test_index]

#now trying loading my data


#need to clean inputs of wrong character
characters_to_remove = ['[', ']', ',']
translation_table = str.maketrans('', '', ''.join(characters_to_remove))


rna_train_data_file = [file for file in input_files if 'rna_train_raw' in file][0]
rna_train_data_file = rna_train_data_file.translate(translation_table).strip() 

protein_train_data_file = [file for file in input_files if 'protein_train_raw' in file][0]
protein_train_data_file = protein_train_data_file.translate(translation_table).strip() 

rna_test_data_file =  [file for file in input_files if 'rna_test_raw' in file][0]
rna_test_data_file = rna_test_data_file.translate(translation_table).strip() 

protein_test_data_file =  [file for file in input_files if 'protein_test_raw' in file][0]
protein_test_data_file = protein_test_data_file.translate(translation_table).strip() 


train_rna = sc.read_h5ad(rna_train_data_file)

#subsample for fewshot
nsamp = int(float(args.ncells))
sampled_obs = train_rna.obs.sample(n=nsamp)
sampled_indices = sampled_obs.index
train_rna_sampled =  train_rna[sampled_indices, :]
train_rna_sampled

test_rna = sc.read_h5ad(rna_test_data_file)

train_protein = sc.read_h5ad(protein_train_data_file)
train_protein_sampled =  train_protein[sampled_indices, :]
test_protein = sc.read_h5ad(protein_test_data_file)

#make sure the current batch size doesn't leave any single-sample batches. 

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
optimal_batch_size_test = find_optimal_batch_size(num_test_samples, args.test_batch_size)
num_train_samples = train_rna_sampled.n_obs  # Number of cells, not features
optimal_batch_size_train = find_optimal_batch_size(num_train_samples, args.batch_size)


if num_test_samples % optimal_batch_size_train == 1:
    optimal_batch_size_train += 1

if num_test_samples % optimal_batch_size_test == 1:
    optimal_batch_size_test += 1

#---  Construct Dataloader ---#
train_kwargs = {'batch_size': args.batch_size}
#test_kwargs = {'batch_size': optimal_batch_size_test}
if use_cuda:
    cuda_kwargs = {'num_workers': 16,
                   'shuffle': False,
                   'prefetch_factor': 2,
                   'pin_memory': True}
    train_kwargs.update(cuda_kwargs)
    #test_kwargs.update(cuda_kwargs)

    my_trainset = fix_SCDataset(train_rna_sampled, train_protein_sampled, args.enc_max_seq_len, args.dec_max_seq_len)



train_sampler = torch.utils.data.distributed.DistributedSampler(my_trainset)
#test_sampler = torch.utils.data.distributed.DistributedSampler(my_testset)

train_loader = torch.utils.data.DataLoader(my_trainset, **train_kwargs, drop_last=False, sampler=train_sampler)
#test_loader = torch.utils.data.DataLoader(my_testset, **test_kwargs, drop_last=False,  sampler=test_sampler)
print("end distributed data")

###############################
#---  Training and Testing ---#
###############################

cell_names = train_protein.obs.index.tolist()

# Create a hash of the cell names
cell_names_str = ''.join(sorted(cell_names))  # Sort for consistency
hash_object = hashlib.md5(cell_names_str.encode())
cell_hash = hash_object.hexdigest()[:24]  # Use first 8 characters

#check for the existence of an already trained model

# Create model filename with the hash
model_filename = f'model_{cell_hash}.pt'
#model_savepath = os.path.join(args.launchdir + "/output/scTranslator/" + model_filename)
model_savepath=model_filename
# Check if model already exists
if os.path.exists(model_filename):
    model = torch.load(model_savepath)
else:
    start_time = time.time()
    for epoch in range(start_epoch+1, args.epochs + 1):
        train_sampler.set_epoch(epoch)
        #test_sampler.set_epoch(epoch)
        torch.cuda.empty_cache()
        
        train_loss, train_ccc = train(args, model, device, train_loader, optimizer, epoch)
        scheduler.step()
        #hash the cell names of the train set, use that as file name
        
    torch.save(model.state_dict(),model_savepath )
 
    
my_testset = fix_SCDataset(test_rna, test_protein, args.enc_max_seq_len, args.dec_max_seq_len)
test_loader = torch.utils.data.DataLoader(my_testset, batch_size=optimal_batch_size_test, drop_last=False)
  
test_loss, test_ccc, test_pcor, y_hat, y, indices = test(model, device, test_loader)


test_cells = test_protein.obs.index
test_feat = test_protein.var['orig_feature']

y_pred =  pd.DataFrame(y_hat,columns=test_feat, index=test_cells)
y_truth = pd.DataFrame(y,columns=test_feat, index=test_cells)
indices = pd.DataFrame(indices)

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
#set rownames and colnames




##############################
#---  Prepare for Storage ---#
##############################

# save results in the first rank
if args.RNA_path == 'dataset/test/dataset1/GSM5008737_RNA_finetune_withcelltype.h5ad':
    dataset_flag = '/seuratv4_16W'
else:
    dataset_flag = '/new_data'
file_path = 'result/test'+dataset_flag
if not os.path.exists(file_path):
    os.makedirs(file_path)
# save args
if rank == 0:
    dict = vars(args)
    filename = open(file_path+'/args'+str(args.repeat)+'.txt','w')
    for k,v in dict.items():
        filename.write(k+':'+str(v))
        filename.write('\n')
    filename.close()

#---  Save the Final Results ---#
# log_path = file_path+'/'+str(rank)+'_rank_log.csv'
# log_all = pd.DataFrame(columns=['train_loss', 'train_ccc', 'test_loss', 'test_ccc'])
# log_all.loc[args.repeat] = np.array([train_loss, train_ccc, test_loss, test_ccc])
# log_all.to_csv(log_path)
# y_pred.to_csv(file_path+'/y_pred.csv')
# y_truth.to_csv(file_path+'/y_truth.csv')
# print('-'*40)
# print('single cell '+str(args.enc_max_seq_len)+' RNA To '+str(args.dec_max_seq_len)+' Protein on dataset'+dataset_flag)
# print('Overall performance on rank_%d in repeat_%d costTime: %.4fs' % (rank, args.repeat, time.time() - start_time))
# print('Training Set: AVG mse %.4f, AVG ccc %.4f' % (np.mean(log_all['train_loss'][:args.repeat]), np.mean(log_all['train_ccc'][:args.repeat])))
# print('Test Set: AVG mse %.4f, AVG ccc %.4f' % (np.mean(log_all['test_loss'][:args.repeat]), np.mean(log_all['test_ccc'][:args.repeat])))



if dobenchmark=='true':
    #get benchmarking prefix
    #make sure looking at filename only, no extended file path
    basename=os.path.basename(rna_train_data_file)
    prefix= basename.replace("_training_data_rna_raw.csv", "") 
    
    log_path = prefix + str(rank)+  '_' + str(args.ncells) + 'cells_fewshot_rank_log.csv'
    #log_path = args.launchdir + '/output/scTranslator/' + prefix + str(rank)+'_finetune_rank_log.csv'
    log_all = pd.DataFrame(columns=['train_loss', 'train_ccc', 'test_loss', 'test_ccc'])
    log_all.loc[args.repeat] = np.array([train_loss, train_ccc, test_loss, test_ccc])
    log_all.to_csv(log_path)
    # 
    #y_pred.to_csv(args.launchdir + "/output/scTranslator/" + prefix + "scTranslator_finetune_prediction.csv")
    #y_truth.to_csv(args.launchdir + "/output/scTranslator/" + prefix + "scTranslator_finetune_truth.csv")
    test_pcor.to_csv(prefix +  '_' + str(args.ncells) + 'cells_fewshot_pearson_corr.csv')
    # test_ccc.to_csv(prefix + '_fewshot_ccc.csv')
    # test_loss.to_csv(prefix + '_fewshot_loss.csv')
    
    y_pred.to_csv( prefix +  '_' + str(args.ncells) + "cells_scTranslator_fewshot_prediction.csv")
    y_truth.to_csv( prefix + '_' + str(args.ncells) +  "cells_scTranslator_fewshot_truth.csv")
    indices.to_csv(prefix +  '_' + str(args.ncells) + "cells_scTranslator_cell_indices.csv")
else:
  
    #log_path = args.launchdir + '/output/scTranslator/' + str(rank)+'_finetune_rank_log.csv'
    log_path =  str(rank)+ '_' + str(args.ncells) +'cells_fewshot_rank_log.csv'
    log_all = pd.DataFrame(columns=['train_loss', 'train_ccc', 'test_loss', 'test_ccc'])
    log_all.loc[args.repeat] = np.array([train_loss, train_ccc, test_loss, test_ccc])
    log_all.to_csv(log_path)
    #y_pred.to_csv(args.launchdir + "/output/scTranslator/scTranslator_finetune_prediction.csv")
    #y_truth.to_csv(args.launchdir + "/output/scTranslator/scTranslator_finetune_truth.csv")
    test_pcor.to_csv(str(args.ncells)  + 'cells_fewshot_pearson_corr.csv')
    # test_ccc.to_csv( 'fewshot_ccc.csv')
    # test_loss.to_csv( 'fewshot_loss.csv')
    
    #pass downstream to pipe
    y_pred.to_csv( str(args.ncells)  + "cells_scTranslator_fewshot_prediction.csv")
    y_truth.to_csv(  str(args.ncells) + "cells_scTranslator_fewshot_truth.csv")
    indices.to_csv(str(args.ncells)  +  "cells_scTranslator_cell_indices.csv")
    
