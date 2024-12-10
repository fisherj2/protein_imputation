
params.scMMT_yaml_dir =  "${moduleDir}/env/scMMT_env.yml"

process scMMT_init_conda{
  label 'small_mem'
  if ( params.do_scMMT) {
    //path to desired conda environment
    conda params.scMMT_yaml_dir 
  }

  input:
   stdin

  output:
  stdout


  script: 
  """
   echo "building scMMT_env"
  """
  
}


//this process defines scipenn method training
process scMMT_train{
  label 'big_mem'
  //path to desired conda environment
  conda params.scMMT_yaml_dir 


  input:

      val file 

  output:

      path "*_prediction.csv"

  script: 
  """
   eval "python  ${moduleDir}/bin/predictADT_scMMT.py --basedir $projectDir --bench ${params.dobenchmark} --files  '${file}'"
  """
}

