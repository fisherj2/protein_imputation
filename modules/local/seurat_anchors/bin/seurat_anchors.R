
#make sure libraries are being loaded from conda env, not any other local R installation
condadir <- Sys.getenv('CONDA_PREFIX')
libpath <- paste0(condadir,'/lib/R/library')
#assign(".lib.loc", libpath, envir = environment(.libPaths))



library(Seurat)
library(Matrix)
library(stringr)
library(data.table)
library(parallel)
set.seed(123)

#Get parameters
args <- commandArgs(trailingOnly = TRUE)
print(args)
basedir<-args[1]
dobenchmark <- args[2]

#clean up file names
args.files <- args[-c(1,2)]

#if remaining args are only length 1, then probably need to split string
if(length(args.files) == 1){
  args.files <- str_split(args.files ,',')[[1]]
}

args.files<-trimws(args.files)
args.files <- str_replace_all(args.files,',|\\[|]','')

#check if output directory exists, create if not
outdir <- paste0(basedir,'/output/Seurat_anchors')

if(!file.exists(outdir)){
  dir.create(outdir, recursive = T)
}



#load in matrices
metadata_file = args.files[grep('metadata_train', args.files)]

rna_train_data_file = args.files[grep('training_data_rna_norm', args.files)]

prot_train_data_file = args.files[grep('training_data_prot_norm', args.files)]

rna_test_data_file = args.files[grep('testing_data_rna_norm', args.files)]


print(args.files)

#build training and testing objects
rna_train_mat <- as.matrix(fread(rna_train_data_file,nThread=detectCores() - 1),rownames=1)
prot_train_mat <- as.matrix(fread(prot_train_data_file,nThread=detectCores() - 1),rownames=1)
rna_test_mat <- as.matrix(fread(rna_test_data_file,nThread=detectCores() - 1),rownames=1)



# Set up training object
train_obj <- CreateSeuratObject(rna_train_mat)

#set normalised data and blank raw counts
train_obj[['RNA']]@data <- train_obj[['RNA']]@counts
train_obj[['RNA']]@counts <- as(Matrix(nrow = 0, ncol = 0, data = 0), "dgCMatrix")
#add protein assay
train_obj[['protein']] <- CreateAssayObject( data=as.matrix(prot_train_mat))

#compute PCA for reference
DefaultAssay(train_obj) <- 'RNA'

#determine feature to use
if(dim(train_obj[['RNA']])[1] < 2000){
  feat.use <- rownames(train_obj[['RNA']])
}else{
  train_obj <- FindVariableFeatures(train_obj)
  feat.use <- train_obj[['RNA']]@var.features
}

train_obj <- ScaleData(train_obj, features = feat.use)
train_obj <- RunPCA(train_obj, features=feat.use)



#set up testing object
test_obj <- CreateSeuratObject(rna_test_mat)

#set normalised data and blank raw counts
test_obj[['RNA']]@data <- test_obj[['RNA']]@counts
test_obj[['RNA']]@counts <- as(Matrix(nrow = 0, ncol = 0, data = 0), "dgCMatrix")




anchors <- FindTransferAnchors(
  features = rownames(train_obj[['RNA']]),
  reference = train_obj,
  reference.reduction = 'pca',
  query = test_obj,
  normalization.method = "LogNormalize",
  dims = 1:50
)


pred <- MapQuery(
  anchorset = anchors,
  query = test_obj,
  reference = train_obj,
  reference.reduction = 'pca',
  refdata = list(
    predicted_protein = "protein"
  )
)


#save predicted expression
pred_mat <- pred@assays$predicted_protein@data


paste0('saving predictions')
if(dobenchmark){
  #get benchmark prefix from files
  prefix <- str_remove(basename(metadata_file),'_metadata.csv')
    #save for later
  write.csv(pred_mat, paste0(basedir, '/output/Seurat_anchors/',prefix,'_seurat_prediction.csv'))
  #pass to pipeline
  write.csv(pred_mat,paste0( prefix,'_seurat_prediction.csv'))
}else{
  #save for later
  write.csv(pred_mat, paste0(basedir, '/output/Seurat_anchors/seurat_prediction.csv'))
  #pass to pipeline
  write.csv(pred_mat, 'seurat_prediction.csv')
}