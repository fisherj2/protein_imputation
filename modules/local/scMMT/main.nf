
//params.scMMT_yaml_dir =  "${moduleDir}/env/scMMT_env.yml"
params.scMMT_yaml_dir = "/shared-workspace/jfisher2/conda/nextflow_envs/scMMT_env-4f2ea4bcf15037255e24c8d62b9c731c"


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
  publishDir "${launchDir}/output/scMMT", mode: 'copy'

  input:

      val file 

  output:

      path "*_prediction.csv"

  script: 
  """
  which python
  eval "python  ${moduleDir}/bin/predictADT_scMMT.py --basedir $projectDir  --bench ${params.dobenchmark} --files  '${file}'"
  """
}

