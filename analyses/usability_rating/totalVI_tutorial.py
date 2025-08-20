
import os
import tempfile

import anndata
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scrublet as scr
import scvi
import seaborn as sns
import torch


pbmc_ref = scvi.data.pbmcs_10x_cite_seq(save_path='/scratch/jfisher/temp/pbmc_10k_protein_v3.h5ad')

pbmc_query = scvi.data.dataset_10x("pbmc_10k_v3", save_path='/scratch/jfisher/temp')

pbmc_query.obs["batch"] = "PBMC 10k (RNA only)"
# put matrix of zeros for protein expression (considered missing)
pro_exp = pbmc_ref.obsm["protein_expression"]
data = np.zeros((pbmc_query.n_obs, pro_exp.shape[1]))
pbmc_query.obsm["protein_expression"] = pd.DataFrame(
    columns=pro_exp.columns, index=pbmc_query.obs_names, data=data
)

#skipped scrublet step

pbmc_query.var["mt"] = pbmc_query.var_names.str.startswith(
    "MT-"
)  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(pbmc_query, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
pbmc_query = pbmc_query[pbmc_query.obs.pct_counts_mt < 15, :].copy()

pbmc_full = anndata.concat([pbmc_ref, pbmc_query])

pbmc_ref = pbmc_full[
    np.logical_or(pbmc_full.obs.batch == "PBMC5k", pbmc_full.obs.batch == "PBMC10k")
].copy()
pbmc_query = pbmc_full[pbmc_full.obs.batch == "PBMC 10k (RNA only)"].copy()


sc.pp.highly_variable_genes(
    pbmc_ref,
    n_top_genes=4000,
    flavor="seurat_v3",
    batch_key="batch",
    subset=True,
)

pbmc_query = pbmc_query[:, pbmc_ref.var_names].copy()

scvi.model.TOTALVI.setup_anndata(
    pbmc_ref, batch_key="batch", protein_expression_obsm_key="protein_expression"
)
totalvi_ref = scvi.model.TOTALVI(pbmc_ref, use_layer_norm="both", use_batch_norm="none")
totalvi_ref.train()

TOTALVI_LATENT_KEY = "X_totalVI"

pbmc_ref.obsm[TOTALVI_LATENT_KEY] = totalvi_ref.get_latent_representation()
sc.pp.neighbors(pbmc_ref, use_rep=TOTALVI_LATENT_KEY)
sc.tl.umap(pbmc_ref, min_dist=0.4)

totalvi_ref_path = os.path.join('/scratch/jfisher/temp', "pbmc_totalvi_ref")
totalvi_ref.save(totalvi_ref_path, overwrite=True)

scvi.model.TOTALVI.prepare_query_anndata(pbmc_query, totalvi_ref_path)
totalvi_query = scvi.model.TOTALVI.load_query_data(
    pbmc_query,
    totalvi_ref_path,
)

totalvi_query.train(200, plan_kwargs={"weight_decay": 0.0})

pbmc_query.obsm[TOTALVI_LATENT_KEY] = totalvi_query.get_latent_representation()
sc.pp.neighbors(pbmc_query, use_rep=TOTALVI_LATENT_KEY)
sc.tl.umap(pbmc_query, min_dist=0.4)

_, imputed_proteins = totalvi_query.get_normalized_expression(
    pbmc_query,
    n_samples=10,
    return_mean=True,
    transform_batch=["PBMC10k", "PBMC5k"],
)