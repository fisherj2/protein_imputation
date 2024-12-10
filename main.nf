#!/usr/bin/env nextflow


// Define DSL2
nextflow.enable.dsl=2



//import modules
include { cTPnet_train; cTPnet_predict; cTPnet_init_conda } from './modules/local/cTPnet'
include { sciPENN_train; sciPENN_init_conda } from './modules/local/scipenn'
include { seurat_anchors; Renv_init_conda} from './modules/local/seurat_anchors'
include { SPECK; SPECK_init_conda } from './modules/local/SPECK'
include { totalvi; totalvi_init_conda } from './modules/local/totalVI'
include { BABEL_train; BABEL_init_conda } from './modules/local/BABEL'
include { scLinear_train; scLinear_init_conda } from './modules/local/scLinear'
include { scMMT_train; scMMT_init_conda } from './modules/local/scMMT'

//benchmarking parameters
cell_int = Channel.from(params.cell_intervals)
prot_int = Channel.from(params.protein_intervals)
comb_int = cell_int.combine(prot_int)



//set default as empty, ovverriden by command line arg
query_file="empty"

// Pipeline proper
// Now we've done the above processes which initiate the required environments, we can specify the processes that actually do the work

//this process takes the input single cell data and formats it into training and test data to be passed downstream
process preprocess_data{
  label 'big_mem'

  //path to desired conda environment
  conda params.Renv_yaml_dir 
  input:

      val inputfile
      val queryfile
      val intervals
      val check

  output:

      path "*.csv", emit: csv_files


  script: 
  """
   Rscript  $projectDir/bin/R/prepare_training_data.R  $projectDir  ${params.dobenchmark} "${intervals}"   ${inputfile}  ${queryfile} 
  """
}


//this process has capacity to convert to h5ad if need for some downstream python analysis. 
process convert_h5ad{
  label 'medium_mem'
  //path to desired conda environment
  conda params.totalVI_yaml_dir 
  
  input:

      val file 

  output:

      path "*.h5ad", emit: csv_files


  script: 
  """
   python  $projectDir/bin/python/convert_h5ad.py  --basedir $projectDir  --files  ${file}
  """
}

//this process defines end-stage evaluation of predictions.
process eval_predictions{
  label 'medium_mem'
  //path to desired conda environment
  conda params.Renv_yaml_dir 
  
  input:

      val file
      val inputfile 
      val queryfile
      val intervals

  output:


  script: 
  
  """

   Rscript $projectDir/bin/R/evaluate_predictions.R  $projectDir ${params.dobenchmark} "${intervals}" "${queryfile}" "${inputfile}"  "${file}"  
  """
}



workflow {
  
  //set up conda env initialization in a serial rather than parallel way. These processes simply cause creation of needed conda envs one at a time (if not already created)
  //This is necessary due to errors caused when multiple conda environemnts are created in parallel
  Renv_init_conda()
  sciPENN_init_conda(Renv_init_conda.out.collect(flat: false))
  cTPnet_init_conda(sciPENN_init_conda.out.collect(flat: false)) 
  SPECK_init_conda(cTPnet_init_conda.out.collect(flat: false))
  totalvi_init_conda(SPECK_init_conda.out.collect(flat: false))
  BABEL_init_conda(totalvi_init_conda.out.collect(flat: false))
  scLinear_init_conda(BABEL_init_conda.out.collect(flat: false))
  scMMT_init_conda(scLinear_init_conda.out.collect(flat: false))



  //convert input to test / training matrices
  preprocess_data(params.input_file, params.query_file, comb_int, scMMT_init_conda.out.collect(flat: false))
  
  //convert_h5ad(preprocess_data.out.csv_files.toList())
  
  //receive test/training files to channel
  input_ch = preprocess_data.out
  
  //for each method check toggle from config file. If true, then run process. If false, emit empty channel

  //compute scipenn predictions
  if ( params.do_scipenn ) {scipenn_out=sciPENN_train(input_ch)} else {scipenn_out=Channel.empty()}
  //compute seurat anchor-based predictions
  if ( params.do_seurat) {seurat_out=seurat_anchors(input_ch)} else {seurat_out=Channel.empty()}
  //compute speck predictions
  if ( params.do_speck ) {speck_out=SPECK(input_ch)} else {speck_out=Channel.empty()}
  //compute scVI/totalVI predictions
  if ( params.do_totalvi ) {totalvi_out=totalvi(input_ch)} else {totalvi_out=Channel.empty()}
  //compute cTPnet predictions
  if ( params.do_cTPnet ) {

    cTPnet_train(input_ch)

    cTPnet_out=cTPnet_predict(cTPnet_train.out.toList(),input_ch)
  }
  else {
    cTPnet_out= Channel.empty()
  }
  // compute BABEL predictions
  if ( params.do_BABEL ) {BABEL_out=BABEL_train(input_ch)} else {BABEL_out=Channel.empty()}
  //compute scLinear predictions
  if ( params.do_scLinear ) {scLinear_out=scLinear_train(input_ch)}else {scLinear_out=Channel.empty()}
  if ( params.do_scMMT ) {scMMT_out=scMMT_train(input_ch)}else {scMMT_out=Channel.empty()}

  ///example one line ifelse for reference
  //if (condition) doThisMethod(); else doThatMethod();

  //collect all outputs from all methods, even the empty ones.
  all_preds = scipenn_out.concat(scMMT_out,scLinear_out, BABEL_out, seurat_out, speck_out, totalvi_out,  cTPnet_out).collect()
  //evaluate the predictions produced
  eval_predictions(all_preds, input_ch.collect(), comb_int, params.query_file)
  
}