message('sys memory')
system('free -m')
options( warn = -1 )

#make sure libraries are being loaded from conda env, not any other local R installation
condadir <- Sys.getenv('CONDA_PREFIX')
libpath <- paste0(condadir,'/lib/R/library')
#assign(".lib.loc", libpath, envir = environment(.libPaths))



library(caret)
library(Seurat)
options(Seurat.object.assay.version = "v3")
library(data.table)
library(parallel)
library(stringr)
set.seed(123)


#Get parameters
args <- commandArgs(TRUE)
basedir<-args[1]
dobenchmark <- as.logical(args[2])
intervals <- args[3]
input_file<-args[4]
queryfile<-args[5]
print(args)

message('loading input file')
print(input_file)

if(endsWith(tolower(input_file),'.rdata')){
  message('loading .Rdata')
  obj <- load(input_file) %>% get()
}else if(endsWith(input_file,'.rds')){
  message('loading .Rds')
  obj <- readRDS(input_file)
}

#ensure assays are version 4 compliant
for(a in names(obj@assays)){
  if(class(obj[[a]]) == 'Assay5'){
    message(paste0('downgrading v5 assay for ',a))
    print(obj[[a]])
    obj[[paste0(a,'v4')]] <- as(object = obj[[a]],Class = 'Assay')
    DefaultAssay(obj) <- paste0(a,'v4')
    obj[[a]] <- NULL
    obj[[a]] <- obj[[paste0(a,'v4')]] 
    DefaultAssay(obj) <- a
    obj[[paste0(a,'v4')]]  <- NULL
    message('downgrade complete')
  }
}


#clean up object components
obj@assays[["RNA"]]@SCTModel.list <- list()

#check if output directory exists, create if not
outdir <- paste0(basedir,'/output/training_files')

if(!file.exists(outdir)){
  dir.create(outdir, recursive = T)
}

outdir <- paste0(basedir,'/output/testing_files')

if(!file.exists(outdir)){
  dir.create(outdir, recursive = T)
}

#prefix for tracking  multiple sets of filtered data
prefix=''

#if query file is specified, need to load that too and do different preparation process 
querygiven <- F
if(queryfile != 'empty'){
  
  querygiven <- T
  #save files based on provided test data
  message('test data object has been provided. Testing on that and ignoring other config partitions.')
  if(endsWith(tolower(queryfile),'.rdata')){
    message('loading .Rdata')
    testobj <- load(queryfile) %>% get()
  }else if(endsWith(queryfile,'.rds')){
    message('loading .Rds')
    testobj <- readRDS(queryfile)
  }
  #ensure assays are version 4 compliant
  for(a in names(testobj@assays)){
    if(class(testobj[[a]]) == 'Assay5'){
      message(paste0('downgrading v5 assay for ',a))
      print(testobj[[a]])

      testobj[[paste0(a,'v4')]] <- as(object = testobj[[a]],Class = 'Assay')
      DefaultAssay(testobj) <- paste0(a,'v4')
      testobj[[a]] <- NULL
      testobj[[a]] <- testobj[[paste0(a,'v4')]] 
      DefaultAssay(testobj) <- a
      testobj[[paste0(a,'v4')]]  <- NULL
      message('downgrade complete')
    }
  }


  #reduce to shared RNA features
  shared_genes <- intersect(rownames(obj[['RNA']]), rownames(testobj[['RNA']]))


  #split normalised RNA data into sets and save
  
  #normalised rna
  message('saving normalised RNA for training')
  training_mat_rna_norm <- obj[['RNA']]@data[shared_genes,]
  fwrite(as.data.frame(training_mat_rna_norm),nThread = detectCores() - 1 , file=paste0(basedir, '/output/training_files/',prefix,'training_data_rna_norm.csv'), row.names=TRUE)
  fwrite(as.data.frame(training_mat_rna_norm),nThread = detectCores() - 1 , file=paste0(prefix,'training_data_rna_norm.csv'), row.names=TRUE)
  rm(training_mat_rna_norm)
  gc()
    message('saving normalised RNA for testing')
  testing_mat_rna_norm <- testobj[['RNA']]@data[shared_genes,]
  fwrite(as.data.frame(testing_mat_rna_norm),nThread = detectCores() - 1 , file=paste0(basedir, '/output/testing_files/',prefix,'testing_data_rna_norm.csv'), row.names=TRUE)
  fwrite(as.data.frame(testing_mat_rna_norm), nThread = detectCores() - 1 ,file=paste0(prefix,'testing_data_rna_norm.csv'), row.names=TRUE)
  rm(testing_mat_rna_norm)
  gc()
  
  
  
  message('saving raw RNA for training')
  training_mat_rna_raw <- obj[['RNA']]@counts[shared_genes,]
  fwrite(as.data.frame(training_mat_rna_raw),nThread = detectCores() - 1 , file=paste0(basedir, '/output/training_files/',prefix,'training_data_rna_raw.csv'), row.names=TRUE)
  fwrite(as.data.frame(training_mat_rna_raw),nThread = detectCores() - 1 , file=paste0(prefix,'training_data_rna_raw.csv'), row.names=TRUE)
  rm(training_mat_rna_raw)
  gc()
  message('saving raw RNA for testing')
  testing_mat_rna_raw <- testobj[['RNA']]@counts[shared_genes,]
  fwrite(as.data.frame(testing_mat_rna_raw),nThread = detectCores() - 1 , file=paste0(basedir, '/output/testing_files/',prefix,'testing_data_rna_raw.csv'), row.names=TRUE)
  fwrite(as.data.frame(testing_mat_rna_raw), nThread = detectCores() - 1 ,file=paste0(prefix,'testing_data_rna_raw.csv'), row.names=TRUE)
  rm(testing_mat_rna_raw)
  gc()


  message('saving normalised protein for training')
  prot_mat_norm <-  obj[['protein']]@data
  training_mat_prot_norm <- prot_mat_norm
  fwrite(as.data.frame(training_mat_prot_norm), nThread = detectCores() - 1 ,file=paste0(basedir, '/output/training_files/',prefix,'training_data_prot_norm.csv'), row.names=TRUE)
  fwrite(as.data.frame(training_mat_prot_norm), nThread = detectCores() - 1 ,file=paste0(prefix,'training_data_prot_norm.csv'), row.names=TRUE)
  rm(training_mat_prot_norm)
  gc()

  message('saving raw protein for training')
  prot_mat_raw <-  obj[['protein']]@counts
  training_mat_prot_raw <- prot_mat_raw
  fwrite(as.data.frame(training_mat_prot_raw), nThread = detectCores() - 1 ,file=paste0(basedir, '/output/training_files/',prefix,'training_data_prot_raw.csv'), row.names=TRUE)
  fwrite(as.data.frame(training_mat_prot_raw), nThread = detectCores() - 1 ,file=paste0(prefix,'training_data_prot_raw.csv'), row.names=TRUE)
  rm(training_mat_prot_raw)
  gc()
  
  #also save metadata
  message('saving metadata')
  meta_df_train <- obj@meta.data
  fwrite(meta_df_train , file=paste0(basedir, '/output/training_files/',prefix,'metadata_train.csv'), row.names=TRUE)
  fwrite(meta_df_train , file=paste0(prefix,'metadata_train.csv'), row.names=TRUE)
  
  meta_df_test <- testobj@meta.data
  fwrite(meta_df_test , file=paste0(basedir, '/output/testing_files/',prefix,'metadata_test.csv'), row.names=TRUE)
  fwrite(meta_df_test , file=paste0(prefix,'metadata_test.csv'), row.names=TRUE)
}else{
  message('no test object provided, partitioning input data')

  #split into training / test set in 80:20 ratio with each celltype evenly represented
  
  ctmeta <- FetchData(obj, 'celltype')
  ctmeta$sample <- rownames(ctmeta)
  
  #check interval specification
  #cleanup unqanted characters
  #intervals <- str_replace_all(intervals, fixed(" "), "")
  intervals <- str_remove_all(intervals,'[(\")(\\)(\r)(\\])(\\[)]')
  intervals <- str_split(intervals,',')[[1]]
  intervals <- trimws(intervals)
  
  #now check if the specified cell interval is numeric or character. If numeric, want a proportion of cells. If character, want a metadata subset.
  cellint <- intervals[1]
  protint <- intervals[2]
  
  print(cellint)
  print(protint)
  
  #cell interval is a string, check further
  negate<-F
  if(!is.na(as.numeric(cellint))){
    cellint <- as.numeric(cellint)
  
    #figure out what proportion of cells / proteins were requested to use
    ncell <- floor(dim(obj)[2] *  cellint)
    
    reduced_ind <- createDataPartition(ctmeta[,1],  p = cellint)
  }else{
  
    if(substr(cellint,1,1) == '_'){
      #we are trying to do leave one out
      negate<-T
      cellint <- substr(cellint,2,nchar(cellint))
    }
    #seek cellint category in object metadata, and set that to be the whole set of used cells
    check <- apply(obj@meta.data,2,FUN=function(x){
     cellint %in% x
    })
  
    
  
  
    if(sum(check) ==0){
      stop('specified cells not found in metadata')
    }else{
      if(sum(check) > 1){
        message('duplicate levels found in metadata. choosing category containing fewest levels')
        cats <- names(which(check))
        meta <- obj@meta.data[,cats]
        levcount <- apply(meta, 2, FUN=function(x){
          length(unique(x))
        })
        cat <- names(which.min(levcount))[1]
      }else{
        #we've found the metadata category with the requested cells. extract
        cat <- names(which(check))
      }
  
      meta <- obj@meta.data[,cat]
      if(negate){
          reduced_ind <- createDataPartition(ctmeta[,1],  p = 1)
          #set test train accordingly
          ct_test <- which(meta == cellint)
          ct_train <- which(meta != cellint)
      }else{
  
        reduced_ind  <- list(first = which(meta == cellint))
      }
  
    }
    
  }
  
  
  
  #check negate to tell whether celltype negation is being done
  if(negate){
    train_cells <- rownames(ctmeta)[reduced_ind[[1]]][ct_train]
    test_cells <- rownames(ctmeta)[reduced_ind[[1]]][ct_test]
  }else{
    trainind <- createDataPartition(ctmeta[reduced_ind[[1]],1],  p = 0.8 )
    train_cells <- rownames(ctmeta)[reduced_ind[[1]]][trainind$Resample1]
    test_cells <- rownames(ctmeta)[reduced_ind[[1]]][which(!(rownames(ctmeta)[reduced_ind[[1]]] %in% c(train_cells)))]
  }
  
  #check if benchmarking. If so, reduce down to representative sample based on input interval
  if(dobenchmark){
    #set prefix for file
    prefix <- paste0(paste0(intervals, collapse='_'),'_')
  }
  
  
  gc()
  
  #split normalised RNA data into sets and save
  
  #normalised rna
  message('saving normalised RNA for training')
  training_mat_rna_norm <- obj[['RNA']]@data[,train_cells]
  fwrite(as.data.frame(training_mat_rna_norm),nThread = detectCores() - 1 , file=paste0(basedir, '/output/training_files/',prefix,'training_data_rna_norm.csv'), row.names=TRUE)
  fwrite(as.data.frame(training_mat_rna_norm),nThread = detectCores() - 1 , file=paste0(prefix,'training_data_rna_norm.csv'), row.names=TRUE)
  rm(training_mat_rna_norm)
  gc()
  message('saving normalised RNA for testing')
  testing_mat_rna_norm <- obj[['RNA']]@data[,test_cells]
  fwrite(as.data.frame(testing_mat_rna_norm),nThread = detectCores() - 1 , file=paste0(basedir, '/output/testing_files/',prefix,'testing_data_rna_norm.csv'), row.names=TRUE)
  fwrite(as.data.frame(testing_mat_rna_norm), nThread = detectCores() - 1 ,file=paste0(prefix,'testing_data_rna_norm.csv'), row.names=TRUE)
  rm(testing_mat_rna_norm)
  gc()
  
  
  message('saving raw RNA for training')
  training_mat_rna_raw <- obj[['RNA']]@counts[,train_cells]
  fwrite(as.data.frame(training_mat_rna_raw),nThread = detectCores() - 1 , file=paste0(basedir, '/output/training_files/',prefix,'training_data_rna_raw.csv'), row.names=TRUE)
  fwrite(as.data.frame(training_mat_rna_raw),nThread = detectCores() - 1 , file=paste0(prefix,'training_data_rna_raw.csv'), row.names=TRUE)
  rm(training_mat_rna_raw)
  gc()
  message('saving raw RNA for testing')
  testing_mat_rna_raw <- obj[['RNA']]@counts[,test_cells]
  fwrite(as.data.frame(testing_mat_rna_raw),nThread = detectCores() - 1 , file=paste0(basedir, '/output/testing_files/',prefix,'testing_data_rna_raw.csv'), row.names=TRUE)
  fwrite(as.data.frame(testing_mat_rna_raw), nThread = detectCores() - 1 ,file=paste0(prefix,'testing_data_rna_raw.csv'), row.names=TRUE)
  rm(testing_mat_rna_raw)
  gc()
  
  message('saving normalised protein for training')
  prot_mat_norm <-  obj[['protein']]@data
  training_mat_prot_norm <- prot_mat_norm[,train_cells]
  
  fwrite(as.data.frame(training_mat_prot_norm), nThread = detectCores() - 1 ,file=paste0(basedir, '/output/training_files/',prefix,'training_data_prot_norm.csv'), row.names=TRUE)
  fwrite(as.data.frame(training_mat_prot_norm), nThread = detectCores() - 1 ,file=paste0(prefix,'training_data_prot_norm.csv'), row.names=TRUE)
  rm(training_mat_prot_norm)
  gc()
  message('saving normalised protein for testing')
  testing_mat_prot_norm <- prot_mat_norm[,test_cells]
  fwrite(as.data.frame(testing_mat_prot_norm),nThread = detectCores() - 1 , file=paste0(basedir, '/output/testing_files/',prefix,'testing_data_prot_norm.csv'), row.names=TRUE)
  fwrite(as.data.frame(testing_mat_prot_norm), nThread = detectCores() - 1 ,file=paste0(prefix,'testing_data_prot_norm.csv'), row.names=TRUE)
  rm(testing_mat_prot_norm)
  gc()
  
  message('saving raw protein for training')
  prot_mat_raw <-  obj[['protein']]@counts
  training_mat_prot_raw <- prot_mat_raw[,train_cells]
  fwrite(as.data.frame(training_mat_prot_raw), nThread = detectCores() - 1 ,file=paste0(basedir, '/output/training_files/',prefix,'training_data_prot_raw.csv'), row.names=TRUE)
  fwrite(as.data.frame(training_mat_prot_raw), nThread = detectCores() - 1 ,file=paste0(prefix,'training_data_prot_raw.csv'), row.names=TRUE)
  rm(training_mat_prot_raw)
  gc()
  message('saving raw protein for testing')
  testing_mat_prot_raw <- prot_mat_raw[,test_cells]
  fwrite(as.data.frame(testing_mat_prot_raw), nThread = detectCores() - 1 ,file=paste0(basedir, '/output/testing_files/',prefix,'testing_data_prot_raw.csv'), row.names=TRUE)
  fwrite(as.data.frame(testing_mat_prot_raw), nThread = detectCores() - 1 ,file=paste0(prefix,'testing_data_prot_raw.csv'), row.names=TRUE)
  rm(testing_mat_prot_raw)
  gc()
  

  meta_df <- obj@meta.data
  message('saving metadata')
  meta_df_train <- meta_df[train_cells,]
  fwrite(meta_df_train , file=paste0(basedir, '/output/training_files/',prefix,'metadata_train.csv'), row.names=TRUE)
  fwrite(meta_df_train , file=paste0(prefix,'metadata_train.csv'), row.names=TRUE)
  
  meta_df_test <-   meta_df_train <- meta_df[test_cells,]
  fwrite(meta_df_test , file=paste0(basedir, '/output/testing_files/',prefix,'metadata_test.csv'), row.names=TRUE)
  fwrite(meta_df_test , file=paste0(prefix,'metadata_test.csv'), row.names=TRUE)

}