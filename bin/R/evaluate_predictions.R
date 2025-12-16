#make sure libraries are being loaded from conda env, not any other local R installation
condadir <- Sys.getenv('CONDA_PREFIX')
libpath <- paste0(condadir,'/lib/R/library')
assign(".lib.loc", libpath, envir = environment(.libPaths))

message('loading libraries')
library(stringr)
library(ggplot2)
library(tidyverse)
library(Hmisc)
library(readr)
#Get parameters
args <- commandArgs(TRUE)
print(args)

basedir<-args[1]
dobenchmark <- as.logical(args[2])
intervals <- args[3]
queryset <- args[4]

args.files <- args[-c(1,2,3,4)]
print(args.files)

#join all remaining files
args.files <- paste0(args.files,collapse=',')

#if remaining args are only length 1, then probably need to split string
if(length(args.files) == 1){
  args.files <- str_split(args.files ,',')[[1]]
}

args.files<-trimws(args.files)
args.files <- str_replace_all(args.files,',|\\[|]','')
print(args.files)

preds <- args.files

#read all the predictions into a list
predList <- list()

#---load in files---

#if benchmarking, need to fetch prefix. filter inputs down to match prefix
if(dobenchmark){
  message('getting intervals')
  intervals <- str_remove_all(intervals,'[(\")(\\)(\r)(\\])(\\[)]')
  intervals <- str_split(intervals,',')[[1]]
  intervals <- trimws(intervals)
  cellint <- intervals[1]
  protint <- intervals[2]

  ind <- grep(paste0('^',cellint), basename(args.files))
  args.files <-args.files[ind]
  message('have reduced arguments to match filter index')
}


print(args.files)

message('loading training and test data')

#record cell and feature names, pulled from training and testing data
ind <- grep('training_data_prot_norm', args.files)[1]
training_data_prot_norm <- read_csv(args.files[ind], col_names = T) %>% as.data.frame()
print(training_data_prot_norm[1:5,1:5])
rownames(training_data_prot_norm) <- training_data_prot_norm[,1]
training_data_prot_norm <- training_data_prot_norm[,-1]
feat <- rownames(training_data_prot_norm)
ind <- grep('testing_data_rna_norm', args.files)[1]
testing_data_rna_norm <- read_csv(args.files[ind], col_names= T) %>% as.data.frame()
rownames(testing_data_rna_norm) <- testing_data_rna_norm[,1]
testing_data_rna_norm<- testing_data_rna_norm[,-1]
cells <- colnames(testing_data_rna_norm)

message('data loaded')

#check for seurat
ind <- grep('seurat', args.files)[1]
if(!is.na(ind) && length(ind) > 0){
  message('loading seurat')
  print(args.files[ind])
  seurat_pred <- read_csv(args.files[ind], col_names=T)
  seurat_pred <- as.data.frame(seurat_pred)
  rownames(seurat_pred) <- seurat_pred[,1]
  seurat_pred <- seurat_pred[,-1]
  colnames(seurat_pred) <- cells
}else{
  seurat_pred <- NULL
}
predList[['Seurat']] <- seurat_pred



# #check for caret file
# ind <- grep('caret', args.files)
# 
# if(length(ind)>0){
#   message('loading caret')
#   caret_preds  <- lapply(ind, FUN=function(i){
#     print(i)
#     df <- as.data.frame(read_csv(args.files[i],col_names = FALSE))
# 
#     #check which way round the matrix is
#     if(is.na(df[1,1])){
#       #first row is colnames and first col is rownames
#       df <- df[-1,-1]
#     }
# 
#     #fix and row/colnames
#     if(dim(df)[1] == length(cells) + 1 || dim(df)[1] == length(feat) + 1 ){
#       rownames(df) <- df[,1]
#       df <- df[,-1]
#     }
# 
#     if(dim(df)[2] == length(cells) + 1 || dim(df)[2] == length(feat) + 1 ){
#       colnames(df) <- df[1,]
#       df <- df[-1,]
#     }
# 
#     if(dim(df)[1] == length(cells) && dim(df)[2]==length(feat)){
#       df<-t(df)
#     }
# 
#     df<-t(apply(df,1,as.numeric))
#     colnames(df) <- cells
#     rownames(df) <- feat
# 
#     return(df)
#     })
# 
#   caret_mods <- lapply(args.files[ind], FUN=function(file){
#     chunk <- str_split(file,'/')[[1]]
#     chunk <- chunk[length(chunk)]
#     chunk <- str_split(chunk,'_')[[1]][2]
#     return(chunk)
#   })
#   names(caret_preds) <- unlist(caret_mods)
#   for(n in names(caret_preds)){
#     predList[[paste0('caret_',n)]] <- caret_preds[[n]]
#   }
# 
# }else{
#   caret_preds <- NULL
#   predList[['caret']] <- caret_preds
# }




ind <- grep('sciPENN', args.files)
if(!is.na(ind) && length(ind) > 0){
  message('loading sciPENN')
  scipenn_pred  <- t(read_csv(args.files[ind], col_names=F))
  rownames(scipenn_pred) <- feat
  colnames(scipenn_pred) <- cells
}else{
  scipenn_pred <- NULL
}
predList[['sciPENN']] <- scipenn_pred

#check for BABEL file

ind <- grep('BABEL', args.files)
if(!is.na(ind) && length(ind) > 0){
  message('loading BABEL')
  BABEL_pred  <- t(read_csv(args.files[ind], col_names=F))
  rownames(BABEL_pred) <- feat
  colnames(BABEL_pred) <- cells
}else{
  BABEL_pred <- NULL
}
predList[['BABEL']] <- BABEL_pred


#check for scMMT file
ind <- grep('scMMT', args.files)
if(!is.na(ind) && length(ind) > 0){
  message('loading scMMT')
  scMMT_pred  <- t(read_csv(args.files[ind], col_names=F))
  rownames(scMMT_pred) <- feat
  colnames(scMMT_pred) <- cells
}else{
  scMMT_pred <- NULL
}
predList[['scMMT']] <- scMMT_pred


#check for scLinear
ind <- grep('sclinear', args.files)
if(!is.na(ind) && length(ind) > 0){
  message('loading scLinear')
  scLinear_pred   <- t(read_csv(args.files[ind], col_names=F))
  colnames(scLinear_pred) <- cells
  rownames(scLinear_pred) <- feat
}else{
  scLinear_pred <- NULL
}
predList[['scLinear']] <- scLinear_pred


#check for totalvi file
ind <- grep('totalVI', args.files)
if(!is.na(ind) && length(ind) > 0){
  message('loading totalVI')
  totalVI_pred  <- t(read_csv(args.files[ind],  col_names=F))
  print(dim(totalVI_pred ))
  print(length(cells))
  print(length(feat))
  colnames(totalVI_pred) <- cells
  rownames(totalVI_pred) <- feat
}else{
  totalVI_pred <- NULL
}
predList[['totalVI']] <- totalVI_pred

#check for SPECK file
ind <- grep('SPECK', args.files)
if(!is.na(ind) && length(ind) > 0){
  message('loading SPECK')

  #speck doesn't use a protein references, only RNA. We need to match the RNA features to our protein measurements in the CITEseq reference. If a reference isn't available, don't evaluate speck
  
  reffile <- paste0(basedir,'/source_data/protein_gene_symbol_map.csv')

  
  if(!file.exists(reffile)){
    SPECK_pred <- NULL
  }else{
    message(paste0( 'loading reference file ',reffile  ) )
    ref <- read.csv(reffile)
    message(paste0('loading SPECK predictions  ',args.files[ind]))
    
    SPECK_pred <- read_csv(args.files[ind])
    SPECK_pred <- as.data.frame(SPECK_pred )
    ind <- match( ref$gene, colnames(SPECK_pred))
    prots <- ref$protein[which(!is.na(ind))]
    
    
    SPECK_pred <- SPECK_pred[,c(1,na.omit(ind))]
    colnames(SPECK_pred) <- prots
    SPECK_pred <- as.matrix(SPECK_pred)
    rownames(SPECK_pred) <- SPECK_pred[,1]
    SPECK_pred <- SPECK_pred[,-1]
    SPECK_pred <- apply(SPECK_pred,1, FUN=function(x){
      as.numeric(x)
    })
    rownames(SPECK_pred) <- prots
    

  }
}else{
  SPECK_pred <- NULL
}

predList[['SPECK']] <- SPECK_pred

#check for  cTP_net
ind <- grep('cTPnet', args.files)
if(!is.na(ind) && length(ind) > 0){
  message('loading cTPnet')
  cTP_net_pred <- read_csv(args.files[ind], col_names = T)
  cTP_net_pred <- as.data.frame(cTP_net_pred)
  rownames(cTP_net_pred) <- cTP_net_pred[,1]
  cTP_net_pred <- cTP_net_pred[,-1]
  rownames(cTP_net_pred) <- feat
  colnames(cTP_net_pred) <- cells
  #fix conversion of periods in feature names
  
}else{
  cTP_net_pred <- NULL
}
predList[['cTP_net']] <- cTP_net_pred

ind <- grep('scTranslator_nofinetune_prediction', args.files)
if(!is.na(ind) && length(ind) > 0){
  message('loading scTranslator_notrain')

  scTranslator_notrain_pred <- read_csv(args.files[ind])
  scTranslator_notrain_pred <- as.data.frame(scTranslator_notrain_pred)
  rownames(scTranslator_notrain_pred) <- scTranslator_notrain_pred [,1]
  scTranslator_notrain_pred <- scTranslator_notrain_pred[,-1]
  scTranslator_notrain_pred <- t(scTranslator_notrain_pred)
  
}else{
  scTranslator_notrain_pred <- NULL
}
predList[['scTranslator_notrain']] <- scTranslator_notrain_pred


ind <- grep('scTranslator_finetune_prediction', args.files)
if(!is.na(ind) && length(ind) > 0){
  message('loading scTranslator_notrain')
  scTranslator_pred <- read_csv(args.files[ind])
  scTranslator_pred <- as.data.frame(scTranslator_pred)
  rownames(scTranslator_pred) <- scTranslator_pred [,1]
  scTranslator_pred <- scTranslator_pred[,-1]
  scTranslator_pred <- t(scTranslator_pred)
  
}else{
  scTranslator_pred <- NULL
}

predList[['scTranslator']] <- scTranslator_pred

ind <- grep('scTranslator_fewshot_prediction', args.files)
if(!is.na(ind) && length(ind) > 0){
  message('loading scTranslator_fewshot')
  scTranslator_pred <- read_csv(args.files[ind])
  scTranslator_pred <- as.data.frame(scTranslator_pred)
  rownames(scTranslator_pred) <- scTranslator_pred [,1]
  scTranslator_pred <- scTranslator_pred[,-1]
  scTranslator_pred <- t(scTranslator_pred)
  
}else{
  scTranslator_pred <- NULL
}

predList[['scTranslator_fewshot']] <- scTranslator_pred


#function for normalizing specific to sctranslator
normalise_scTranslator <- function(x){
  MIN = min(x)
  MAX = max(x)
  newx = 1e-8 +  (x-MIN) / (MAX - MIN )  * (1 - 1e-8)
}

corrList <- list()
RMSEList <- list()

#load true protein values if they exist. If no testing protein data exists, we are using an RNA only query and can't check for RMSE / prediction correlation
#else if a separate training and query set were provided at command line, don't do this step either. 
ind1 <- grep('testing_data_prot_raw', args.files)[1]
ind2 <- grep('testing_data_prot_norm', args.files)[1]

if(  length(ind1) > 0 && length(ind2) > 0 && queryset == 'empty' && !is.na(ind1) && !is.na(ind2)){
  
  #we are training and testing on the same dataset, can compute correlation and RMSE with self
  testing_data_prot_raw <- read.csv(args.files[ind1], row.names = 1)
  testing_data_prot_norm <- read.csv(args.files[ind2], row.names = 1)

  message('computing metrics')
  #compute metrics between the predictions and the true protein

  for(n in names(predList)){
    print(n)
    pred <- predList[[n]]

    if(grepl('scTranslator',n)){

      #need to compute w.r.t. the special normalisation used by scTranslator
      testing_data_prot_norm_minmax <-apply(testing_data_prot_raw,2,FUN=normalise_scTranslator) 
        
      corVal <- cor(t(pred), t(testing_data_prot_norm_minmax[rownames(pred),colnames(pred)] ), method = 'spearman')
    }else{
      corVal <- cor(t(pred), t(testing_data_prot_norm[rownames(pred),] ), method = 'spearman')
    }
    
    corVal <- diag(corVal)
    names(corVal) <- rownames(pred)
    corrList[[n]] <- corVal

    #rmse
    if(grepl('scTranslator',n)){
      diff <- pred - testing_data_prot_norm_minmax[rownames(pred),]
    }else{
      diff <- pred - testing_data_prot_norm[rownames(pred),]
    }
    
    RMSEval <- apply(diff, 1, FUN=function(x){sqrt(mean(x)^2)})
    RMSEList[[n]] <- RMSEval
  }

}


message('saving outcomes')
if(dobenchmark){
  ind <- grep('testing_data_prot_norm', args.files)
  filename <- args.files[ind]

  chunks <- str_split(basename(filename),'_')[[1]][c(1,2)]
  prefix <- paste0(chunks,collapse='_')
  save(predList,corrList, RMSEList,  file= paste0(prefix,'_bechmark_predictions.Rdata'))
}else{
  save(predList,corrList, RMSEList,  file= 'predictions.Rdata')
}




