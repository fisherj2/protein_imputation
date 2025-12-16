
//params.cTPnet_yaml_dir =  "${moduleDir}/env/cTPnet_exported_env.yml"
params.cTPnet_yaml_dir = "/shared-workspace/jfisher2/conda/nextflow_envs/cTPnet_exported_env-a73e060ea1874b0e106be82269e76db0"


process cTPnet_init_conda{
  label 'small_mem'
  if ( params.do_cTPnet) {
    //path to desired conda environment
    conda params.cTPnet_yaml_dir
  }


  input:
   stdin

  output:
   stdout


  script: 
  """
   echo "building ctpnet_env"
  """
  
}


//this process defines cTPnet model training
process cTPnet_train{
  label 'big_mem'
  //path to desired conda environment
  conda params.cTPnet_yaml_dir 
  publishDir "${launchDir}/output/cTPnet", mode: 'copy'
  
  input:

      val file 

  output:

      path '*final_model_rep*'    

  script: 
  """
   echo 'cTPnet training'
   eval "python -s ${moduleDir}/bin/cTP_net_train.py  --basedir $projectDir --bench ${params.dobenchmark} --files  '${file}'"
  """
}

// this process takes input from the above (python), and uses it to predict on test data (R)
process cTPnet_predict{
  label 'big_mem'
  //path to desired conda environment
  conda params.cTPnet_yaml_dir 
  
  input:
      val model
      val file 

  output:

      path "*_prediction.csv"

  script: 
  """
   echo 'running cTPnet prediction function'
   eval "Rscript  ${moduleDir}/bin/cTPnet_predict.R   $projectDir  ${params.dobenchmark} '${moduleDir}'  '${model}'  '${file}' "
  """
}