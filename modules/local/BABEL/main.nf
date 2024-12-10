
params.BABEL_yaml_dir =  "${moduleDir}/env/BABEL_env.yml"

process BABEL_init_conda{
  label 'small_mem'
  if ( params.do_BABEL) {
    //path to desired conda environment
    conda params.BABEL_yaml_dir 
  }

  input:
   stdin

  output:
    stdout


  script: 
  """
   echo "building BABEL_env"
  """
  
}


//this process defines scipenn method training
process BABEL_train{
  label 'big_mem'
  //path to desired conda environment
  conda params.BABEL_yaml_dir 
  

  input:

      val file 

  output:

      path "*_prediction.csv"

  script: 
  """
   eval "python  ${moduleDir}/bin/predictADT_BABEL.py --basedir $projectDir  --bench ${params.dobenchmark} --files  '${file}' "
  """
}

