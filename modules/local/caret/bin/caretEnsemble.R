
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
dobenchmark <- args[2]

#identify models provided
model_params <- args[3]

model_params <- str_split(model_params,',')[[1]]

model_params <- unlist(lapply(model_params, FUN=function(x){
  str_replace_all(x,' |,|\\[|]','')
}))

modelType <- model_params[1]
responsepath <- model_params[2]

print('response')
print(responsepath)
#everything else is rest of files
args.files <- args[-c(1,2,3)]
args.files <- str_replace_all(args.files,' |,|\\[|]','')


#check if output directory exists, create if not
outdir <- paste0(basedir,'/output/caret')

if(!file.exists(outdir)){
  dir.create(outdir, recursive = T)
}

# #prep testing full data
# load("/fcrbiouatappn01/resbioinfo/data/MiroBio/shared/home/jfisher/analyses/protein_prediction/formal_analysis/source_data/Hao_et_al_Seurat_obj.Rdata")
# rna_train_mat <- as.matrix(reference[['RNA']]@data[,sample(1:dim(reference)[2], round(dim(reference)[2] * 0.8))])
# prot_train_mat <- as.matrix(reference[['protein']]@data[,colnames(rna_train_mat)])
#
# basedir = '/fcrbiouatappn01/resbioinfo/data/MiroBio/shared/home/jfisher/analyses/protein_prediction/formal_analysis'
# metadata_file =  '/scratch/jfisher2/protein_prediction/test_output/training_files/metadata.csv'
# rna_train_data_file = '/scratch/jfisher2/protein_prediction/test_output/training_files/training_data_rna_norm.csv'
# prot_train_data_file = '/scratch/jfisher2/protein_prediction/test_output/training_files/training_data_prot_norm.csv'
# rna_test_data_file =  '/scratch/jfisher2/protein_prediction/test_output/testing_files/testing_data_rna_norm.csv'
# prot_test_data_file =  '/scratch/jfisher2/protein_prediction/test_output/testing_files/testing_data_prot_norm.csv'


#load in matrices
metadata_file = args.files[grep('metadata', args.files)]
rna_train_data_file = args.files[grep('training_data_rna_norm', args.files)]
prot_train_data_file = responsepath
rna_test_data_file = args.files[grep('testing_data_rna_norm', args.files)]
prot_test_data_file = args.files[grep('testing_data_prot_norm', args.files)]

metadata <- as.matrix(fread(metadata_file,nThread=detectCores() - 1),rownames=1)
rna_train_mat <- as.matrix(fread(rna_train_data_file,nThread=detectCores() - 1),rownames=1)
prot_train_mat <- as.matrix(fread(prot_train_data_file,nThread=detectCores() - 1),rownames=1)
rna_test_mat <- as.matrix(fread(rna_test_data_file,nThread=detectCores() - 1),rownames=1)
prot_test_mat <- as.matrix(fread(prot_test_data_file,nThread=detectCores() - 1),rownames=1)


test_dat <- t(rna_test_mat)
rm(rna_test_mat)
gc()

train_dat <- t(rna_train_mat)
train_dat <- apply(train_dat,2,as.numeric)

if(dim(train_dat)[2] > 500){
  print('reducing to 500 most variable features')
  #clean zero variance features
  #nzv <- nearZeroVar(train_dat, saveMetrics= TRUE,allowParallel=TRUE, foreach=T)
  #train_dat <- train_dat[, rownames(nzv)[-which(nzv$zeroVar)]]

  #take top 200 variable features
  varvals<-FindVariableFeatures(t(train_dat))
  varvals <- varvals[order(varvals$variance.standardized, decreasing = T),]
  keepfeat<-rownames(varvals[1:min(500, dim(varvals)[1]),])
  train_dat<-train_dat[,keepfeat]
  test_dat <- test_dat[,keepfeat]
}

rm(rna_train_mat)

gc()


#set number of cores so each one has enough ram to hold data and run computation on it

#work out sensible data size to cpu ratio
datasize <- print(object.size(train_dat) + object.size(prot_train_mat) + object.size(test_dat) + object.size(prot_test_mat),units="Mb", standard='SI')
compsize <- datasize *3
neededsize <- datasize + compsize
freeram <- benchmarkme::get_ram()

core_num <- min( detectCores() - 1 , as.integer(freeram / neededsize))

print(paste0('running caret with ', core_num , ' cores'))
cl <- makePSOCKcluster(core_num)
registerDoParallel(cl)



#set parameters
xgbTree_grid <- expand.grid(
  nrounds = 100,
  max_depth = 6,
  eta = 0.3,
  gamma = 0,
  colsample_bytree = 1,
  min_child_weight = 1,
  subsample = 1
)

rpart_grid <- expand.grid(
  cp=0.01
)

ranger_grid <- expand.grid(
  mtry = sqrt(ncol(train_dat)),
  splitrule = 'variance',
  min.node.size  = 5
)

tunelist <- list(
                  rpart=caretModelSpec(method="rpart", tuneGrid = rpart_grid),
                  glm=caretModelSpec(method="glm"),
                  lm=caretModelSpec(method="lm"),
                  xgbTree=caretModelSpec(method="xgbTree", tuneGrid = xgbTree_grid),
                  rf=caretModelSpec(method="ranger", tuneGrid=ranger_grid)
                )

#adjust to match requested inputs
tunelist <- tunelist[modelType]


#iterate over each protein
prediction_list <- list()


prot <- rownames(prot_train_mat)[1]
tic()
message(paste0('fitting caret models for ', prot))
message(Sys.time())
protein <- data.frame(protein = prot_train_mat[prot,])
colnames(protein) <-'protein'


my_control <- trainControl("cv",
                           number=5,
                           savePredictions = 'final',
                           verboseIter = TRUE,
                           allowParallel = TRUE
                           )

print('starting model training')
model <- caretList(
            x = train_dat,
            y= protein[,1],
            trControl=my_control,
            #methodList=modelTypes,
            tuneList=tunelist)



#save models to file
model_pred <- lapply(model, predict, newdata=test_dat[,colnames(train_dat)])
toc()




protname <- str_replace_all(prot,'/','-')
if(dobenchmark){
  #get benchmark prefix from files
  chunks <- str_split(rna_train_data_file,'_')[[1]][c(1,2)]
  prefix <- paste0(chunks,collapse='_')
  write.csv(model_pred , file=paste0(basedir,'/output/caret/',prefix,'_caret_', modelType,'_',protname,'_prediction.csv'))
  write.csv(model_pred , file=paste0( prefix,'_caret_', modelType,'_',protname,'_prediction.csv'))
}else{

  write.csv(model_pred , file=paste0(basedir,'/output/caret/','caret_', modelType,'_',protname,'_prediction.csv'))
  write.csv(model_pred , file=paste0( 'caret_', modelType,'_',protname,'_prediction.csv'))
}


message(paste0('have written predictions for ',modelType,  ' ', prot))


message('caret iteration finished')