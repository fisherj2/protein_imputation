
params.sciPENN_yaml_dir =  "${moduleDir}/env/sciPENN_env.yml"

process sciPENN_init_conda{
  label 'small_mem'
  if ( params.do_scipenn) {
    //path to desired conda environment
    conda params.sciPENN_yaml_dir 
  }


  input:
   stdin

  output:
    stdout


  script: 
  """
  if [[ "${params.do_scipenn}" == "true" ]]
  then
    echo "building scipenn_env"

    #check if scipenn installed. if not run installation script. Need to do it like this because dependencies are specified wrong.
    source ${moduleDir}/bin/check_scipenn_env.sh
  fi
  """
  
}


//this process defines scipenn method training
process sciPENN_train{
  label 'big_mem'
  //path to desired conda environment
  conda params.sciPENN_yaml_dir 
  

  input:

      val file 

  output:

      path "*_prediction.csv"

  script: 
  """

    eval "python  ${moduleDir}/bin/scipenn_train.py  --basedir $projectDir --bench ${params.dobenchmark} --files  '${file}'"

  """
}

