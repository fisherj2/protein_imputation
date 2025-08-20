

params.totalVI_yaml_dir= "${moduleDir}/env/scVI_exported_env.yml"



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
  
  input:

      val file 

  output:

      path "*_prediction.csv"

  script: 
  """
   eval "python  ${moduleDir}/bin/scVI_train.py  --basedir $projectDir --launchdir $launchDir  --bench ${params.dobenchmark} --files  '${file}'"
  """
}
