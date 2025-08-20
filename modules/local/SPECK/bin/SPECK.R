#make sure libraries are being loaded from conda env, not any other local R installation
condadir <- Sys.getenv('CONDA_PREFIX')
libpath <- paste0(condadir,'/lib/R/library')
assign(".lib.loc", libpath, envir = environment(.libPaths))



# these dependencies should be installed by the conda .yml but it wasn't working properly so put these checks here as a backup to ensure correct versions are used. 
library(remotes)
if(packageVersion("Matrix") < "1.6.3") {
  install_version(
    'Matrix',
    version = '1.6.3',
    repos='http://cran.us.r-project.org')
}

if(packageVersion("rlang") < "1.1.0") {
  install_version(
    'rlang',
    version = '1.1.0',
    repos='http://cran.us.r-project.org')
}


if(packageVersion("Rcpp") < "1.0.8.3") {
  install_version(
    'Rcpp',
    version = '1.0.8.3',
    repos='http://cran.us.r-project.org')
}

if(packageVersion("generics") < "0.1.3") {
  install_version(
    'generics',
    version = '0.1.3',
    repos='http://cran.us.r-project.org')
}



#check for package and install if missing
library(Matrix)
library(rlang)
library(generics)
library(stringr)


#SPECK only available on CRAN, and depends on Seurat 5 also only on CRAN at time of writing. We install here as cannot install via .yml

if (!require("SPECK")) install.packages("SPECK",  repos='http://cran.us.r-project.org', dependencies=T)
library(SPECK)
library(Seurat)

#Get parameters
args <- commandArgs(trailingOnly = TRUE)
print(args)
basedir<-args[1]
launchdir <- args[2]
dobenchmark <- args[3]

#clean up file names
args.files <- args[-c(1,2,3)]

#if remaining args are only length 1, then probably need to split string
if(length(args.files) == 1){
  args.files <- str_split(args.files ,',')[[1]]
}

args.files<-trimws(args.files)
args.files <- str_replace_all(args.files,',|\\[|]','')


#check if output directory exists, create if not
outdir <- paste0(launchdir,'/output/SPECK')

if(!file.exists(outdir)){
  dir.create(outdir, recursive = T)
}


#read in arguments
metadata_file = args.files[grep('metadata_train', args.files)]
rna_test_data_file = args.files[grep('testing_data_rna_raw', args.files)]


#build training and testing objects

rna_test_mat <- read.csv(rna_test_data_file, row.names = 1)



speck.full <- speck(counts.matrix = t(rna_test_mat), rank.range.end = min(100, dim(rna_test_mat)[2]),
                    min.consec.diff = 0.01, rep.consec.diff = 2,
                    manual.rank = NULL, max.num.clusters = 4,
                    seed.rsvd = 1, seed.ckmeans = 2)
speck.rank <- speck.full$rrr.rank


speck.output<-speck.full$thresholded.mat


if(dobenchmark){
  #get benchmark prefix from files
  prefix <- str_remove(basename(metadata_file),'_metadata.csv')
    #save for later
  write.csv(speck.output, paste0(launchdir, '/output/SPECK/',prefix,'_SPECK_prediction.csv'))
  #pass to pipeline
  write.csv(speck.output,paste0( prefix,'_SPECK_prediction.csv'))
}else{
  #save for later
  write.csv(speck.output, paste0(launchdir, '/output/SPECK/SPECK_prediction.csv'))
  #pass to pipeline
  write.csv(speck.output, 'SPECK_prediction.csv')
}