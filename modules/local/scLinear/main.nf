
params.scLinear_yaml_dir =  "${moduleDir}/env/scLinear_env.yml"

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


  input:

      val file 

  output:

      path "*_prediction.csv"

  script: 
  """
   eval "python  ${moduleDir}/bin/predictADT_scLinear.py --basedir $projectDir --launchdir $launchDir  --bench ${params.dobenchmark} --files  '${file}'"
  """
}

