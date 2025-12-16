

//params.totalVI_yaml_dir= "${moduleDir}/env/scVI_exported_env.yml"

params.totalVI_yaml_dir="/shared-workspace/jfisher2/conda/nextflow_envs/scVI_exported_env-9f010de8e1e51acc476a048ea18e4355"


process totalvi_init_conda{
  label 'small_mem'
  if ( params.do_totalvi) {
    //path to desired conda environment
    conda params.totalVI_yaml_dir
  }


  input:
   stdin

  output:
    stdout


  script: 
  """
   echo "building totalvi_env"

  """
  
}


process totalvi{
  label 'medium_mem'
  //path to desired conda environment
  conda params.totalVI_yaml_dir 
  publishDir "${launchDir}/output/totalVI", mode: 'copy'
  
  input:

      val file 

  output:

      path "*_prediction.csv"
      path "*.pt"

  script: 
  """
   eval "python  ${moduleDir}/bin/scVI_train.py  --basedir $projectDir   --bench ${params.dobenchmark} --files  '${file}'"
  """
}
