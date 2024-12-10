
#make sure libraries are being loaded from conda env, not any other local R installation
condadir <- Sys.getenv('CONDA_PREFIX')
libpath <- paste0(condadir,'/lib/R/library')
#assign(".lib.loc", libpath, envir = environment(.libPaths))


library(caret)
library(rpart)
library(caretEnsemble)
library(caTools)
library(doParallel)
library(stringr)
library(tictoc)
library(Seurat)
library(data.table)
library(parallel)
set.seed(123)

#Get parameters
args <- commandArgs(TRUE)
print(args)

basedir<-args[1]

#clean up file names
args.files <- args[-c(1)]
args.files <- str_replace_all(args.files,' |,|\\[|]','')


#check if output directory exists, create if not
outdir <- paste0(basedir,'/output/caret')

if(!file.exists(outdir)){
  dir.create(outdir, recursive = T)
}

#load in training output
prot_train_data_file = args.files[grep('training_data_prot_norm', args.files)]


prot_train_mat <- as.matrix(fread(prot_train_data_file,nThread=detectCores() - 1),rownames=1)

#partition and save to split across multiple processes

for(prot in rownames(prot_train_mat)){
    dat <- prot_train_mat[prot,,drop=F]

    #fix name
    protname <- prot
    protname <-   str_replace_all(prot, '/','-')
    write.csv(dat, file=paste0( 'caret_', protname,'_training_y.csv'), row.names=TRUE)
}
