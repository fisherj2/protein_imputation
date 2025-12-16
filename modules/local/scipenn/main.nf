
//params.sciPENN_yaml_dir =  "${moduleDir}/env/sciPENN_exported_env.yml"
params.sciPENN_yaml_dir ="/shared-workspace/jfisher2/conda/nextflow_envs/sciPENN_exported_env-e8006b678a02d63d451b1782744f3c68"


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
    source ${moduleDir}/bin/check_scipenn_env.sh
  fi
  """
  
}


//this process defines scipenn method training
process sciPENN_train{
  label 'big_mem'
  //path to desired conda environment
  conda params.sciPENN_yaml_dir 
   publishDir "${launchDir}/output/sciPENN", mode: 'copy'  //

  input:

      val file 
      val check_scipenn_init

  output:

      path "*_prediction.csv"

  script: 
  """

    eval "python  ${moduleDir}/bin/scipenn_train.py  --basedir $projectDir --bench ${params.dobenchmark} --files  '${file}'"

  """
}

