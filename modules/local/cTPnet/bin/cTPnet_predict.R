
# #needed software is only available on github, and some dependencies are only available on CRAN. So we install here. 
# library(devtools)
# if (!require("cTPnet")) devtools::install_github("zhouzilu/cTPnet")
# 
# 
# library(cTPnet)
library(data.table)
library(parallel)

message('starting cTPnet prediction')
#make sure libraries are being loaded from conda env, not any other local R installation
condadir <- Sys.getenv('CONDA_PREFIX')
libpath <- paste0(condadir,'/lib/R/library')
#assign(".lib.loc", libpath, envir = environment(.libPaths))

#set
reticulate::use_condaenv(condadir, required = TRUE)

library(Seurat)
library(reticulate)
library(stringr)

#Get parameters
args <- commandArgs(trailingOnly = TRUE)
basedir<-args[1]
dobenchmark <- args[2]
scriptdir <- args[3]



#clean up file names
models <- args[4]


args.files <- args[-c(1,2,3,4)]



print(args.files)
#if remaining args are only length 1, then probably need to split string
if(length(args.files) == 1){
  message('splitting input')
  args.files <- str_split(args.files ,',')[[1]]
}

if(length(models) == 1){
  message('splitting input')
  models <- str_split(models ,',')[[1]]
}


args.files<-trimws(args.files)
args.files <- str_replace_all(args.files,',|\\[|]','')

models <- trimws(models)
models<-str_replace_all(models,',|\\[|]','')

#remove double slashes
print(models)

print(args.files)

metadata_file = args.files[grep('metadata_train', args.files)]


if(dobenchmark){
  #we need to pay attention to the prefix, to match the right model to the input files
  prefix <- str_remove(basename(metadata_file),'_metadata_train.csv')
  ind <- intersect(grep('final_model', models), grep(paste0('^',prefix), basename(models)))
  modeldir <- models[ind]
}else{
  #just use the first model file we find
  ind <- grep('final_model', models)
  modeldir <- models[ind]
}

print(modeldir)
#check if output directory exists, create if not
outdir <- paste0(basedir,'/output/cTPnet')

if(!file.exists(outdir)){
  dir.create(outdir, recursive = T)
}

message('fetching files')
#fetch input data paths
#load in matrices

rna_train_data_file = args.files[grep('training_data_rna_norm', args.files)]

prot_train_data_file = args.files[grep('training_data_prot_norm', args.files)]

rna_test_data_file = args.files[grep('testing_data_rna_norm', args.files)]



# 
# #temp paths here for testing !!
# modeldir = '/scratch/jfisher2/protein_prediction/output/cTPnet/final_model_rep0_ep65'
# basedir = '/scratch/jfisher2/protein_prediction'
# metadata_file =  '/scratch/jfisher2/protein_prediction/output/training_files/metadata.csv'
# rna_train_data_file = '/scratch/jfisher2/protein_prediction/output/training_files/training_data_rna_raw.csv'
# prot_train_data_file = '/scratch/jfisher2/protein_prediction/output/training_files/training_data_prot_raw.csv'
# rna_test_data_file =  '/scratch/jfisher2/protein_prediction/output/testing_files/testing_data_rna_raw.csv'



#check and records the dimensions of the training data and included protein features
metadata <- read.csv(metadata_file, row.names = 1)
rna_train_mat <- as.matrix(fread(rna_train_data_file,nThread=detectCores() - 1),rownames=1)
prot_train_mat <- as.matrix(fread(prot_train_data_file,nThread=detectCores() - 1),rownames=1)
rna_test_mat <- as.matrix(fread(rna_test_data_file,nThread=detectCores() - 1),rownames=1)

prot_feat <- rownames(prot_train_mat)
rna_data_dim <- dim(rna_train_mat)[1]
rna_feat <- rownames(rna_train_mat)

#clean problematic characters
prot_feat <- str_remove_all(prot_feat,'\\.')
rna_feat <- str_remove_all(rna_feat,'\\.')

#---- PREDICTION FROM MODEL---

#we have modified some functions from the cTPnet github repo, and stored them under /bin
# we source these now to apply prediction on a cTPnet object via R functions which call reticulate
print('loading python prediction function')
source_python(paste0(scriptdir ,'/bin/cTP_net_predict.py'))

source(paste0(scriptdir ,'/bin/cTPnet_postprocess.R'))
source(paste0(scriptdir ,'/bin/cTPnet_preprocess.R'))
source(paste0(scriptdir ,'/bin/modified_cTPnet.R'))


#model_file_path="/export/home/jfisher22/jfisher2/analyses/protein_prediction/source_data/cTPnet_weight_24"

data_type='Seurat3'

print('running prediction')
test_obj <- CreateSeuratObject(rna_test_mat, meta.data=metadata[colnames(rna_test_mat),])
test_obj <-  modified_cTPnet(test_obj,data_type,modeldir, dprotein='custom',  protein_names=prot_feat, rna_names=rna_feat,rna_dim=rna_data_dim)


#save prediction
#save predicted expression
pred_mat <- test_obj@assays$cTPnet@data



paste0('saving predictions')
if(dobenchmark){
  #get benchmark prefix from files
  chunks <- str_split(basename(rna_test_data_file),'_')[[1]][c(1,2)]
  prefix <- paste0(chunks,collapse='_')
    #save for later
  write.csv(pred_mat, paste0(basedir, '/output/cTPnet/',prefix,'_cTPnet_prediction.csv'))
  #pass to pipeline
  write.csv(pred_mat,paste0( prefix,'_cTPnet_prediction.csv'))
}else{
  #save for later
  write.csv(pred_mat, paste0(basedir, '/output/cTPnet/cTPnet_prediction.csv'))
  #pass to pipeline
  write.csv(pred_mat, 'cTPnet_prediction.csv')
}

