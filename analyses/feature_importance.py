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


pipe = ADTPredictor(do_log1p=True)
pipe.load(path="/home/jfisher2/workspace/analyses/protein_prediction_publication/benchmarking/method_predictions/1_1_sclinear_model")


#load in test data
adata_rna_test = ad.read_csv("/home/jfisher2/workspace/analyses/protein_prediction_publication/benchmarking/testing_files/1_1_testing_data_rna_raw.csv", delimiter=',', first_column_names = True).transpose()

impvals = pipe.feature_importance( gex_test= adata_rna_test.X)

np.savetxt("/scratch/jfisher2/protein_prediction/output/scLinear/1_1_model_feature_importance.csv" ,impvals,delimiter=",")
