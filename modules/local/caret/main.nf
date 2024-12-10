

params.caret_yaml_dir =  "${moduleDir}/env/caret_env.yml"



process caret_init_conda{
  label 'small_mem'
  if ( params.do_caret) {
    //path to desired conda environment
    conda params.caret_yaml_dir
  }
  input:
   stdin

  output:
    stdout
  script: 
  """

   echo "building caret_env"

  """
  
}


process caret_split_output{

  label 'medium_mem'
  conda params.caret_yaml_dir
  
  input:

      val file 

  output:

      path "*_training_y.csv" 

  script: 
  """
   echo 'splitting caret training data'
   Rscript ${moduleDir}/bin/splitTraining.R  $projectDir  ${file}
  """
  
}




process caret_train{
  label 'medium_mem'
  conda params.caret_yaml_dir
  
  
  input:
      val training_output
      val file 

  output:

      path "*_prediction.csv" 

  script: 
  """
  echo 'beginning caret training'
   Rscript ${moduleDir}/bin/caretEnsemble.R  $projectDir  ${params.dobenchmark}  "${training_output}" ${file}
  """
}


