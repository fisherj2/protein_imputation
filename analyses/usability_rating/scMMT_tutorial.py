import numpy as np
import pandas as pd
import pickle
from time import time
from scipy.stats import spearmanr, gamma, poisson
from anndata import AnnData, read_h5ad
import scanpy as sc
from scanpy import read
import torch
from torch.utils.data import DataLoader, TensorDataset
from torch import tensor
from torch.cuda import is_available
from scMMT.scMMT_API import scMMT_API
from sklearn.metrics import f1_score, accuracy_score
import warnings
warnings.filterwarnings("ignore")

seed = 5
torch.manual_seed(seed)
torch.cuda.manual_seed(seed) 
torch.cuda.manual_seed_all(seed)  # if you are using multi-GPU
np.random.seed(seed)

# Read in Raw Data
adata_gene = sc.read("/scratch/jfisher2/temp/pbmc_gene.h5ad")
adata_protein = sc.read("/scratch/jfisher2/temp/pbmc_protein.h5ad")

# This is the protein processing process, which can be switched to any processing method
sc.pp.normalize_total(adata_protein)
sc.pp.log1p(adata_protein)
patients = np.unique(adata_protein.obs['donor'].values)
for patient in patients:
    print(patient)
    indices = [x == patient for x in adata_protein.obs['donor']]
    sub_adata = adata_protein[indices]
    sc.pp.scale(sub_adata)
    adata_protein[indices] = sub_adata.X


# Create training and testing dataset
train_bool = [x in ['P1', 'P3', 'P4', 'P7'] for x in adata_protein.obs['donor']]
adata_gene_train = adata_gene[train_bool].copy()
adata_protein_train = adata_protein[train_bool].copy()
adata_gene_test = adata_gene[np.invert(train_bool)].copy()
adata_protein_test = adata_protein[np.invert(train_bool)].copy()


scMMT = scMMT_API(    gene_trainsets = [adata_gene_train], protein_trainsets = [adata_protein_train], gene_test = adata_gene_test, 
                      train_batchkeys = ['donor'], test_batchkey = 'donor',
                      log_normalize = True,            # Is scRNA seq standardized for log
                      type_key = 'celltype.l3',        # Keywords representing cell types (in protein dataset)
                      data_dir="/scratch/jfisher2/temp/preprocess_data_l3.pkl",  # Save path for processed data
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