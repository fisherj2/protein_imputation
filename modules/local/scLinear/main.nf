
//params.scLinear_yaml_dir =  "${moduleDir}/env/scLinear_env.yml"
params.scLinear_yaml_dir = "/shared-workspace/jfisher2/conda/nextflow_envs/scLinear_env-84a3f3e706ec1b57cc8a3579527980fb"


process scLinear_init_conda{
  label 'small_mem'
  if ( params.do_scLinear) {
    //path to desired conda environment
    conda params.scLinear_yaml_dir 
  }

  input:
   stdin

  output:
  stdout


  script: 
  """
   echo "building scLinear_env"
  """
  
}


//this process defines scipenn method training
process scLinear_train{
  label 'big_mem'
  //path to desired conda environment
  conda params.scLinear_yaml_dir 
  publishDir "${launchDir}/output/scLinear", mode: 'copy' 


  input:

      val file 

  output:

      path "*_prediction.csv"

  script: 
  """
   eval "python  ${moduleDir}/bin/predictADT_scLinear.py --basedir $projectDir  --bench ${params.dobenchmark} --files  '${file}'"
  """
}

