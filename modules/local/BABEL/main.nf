
//params.BABEL_yaml_dir =  "${moduleDir}/env/BABEL_exported_env.yml"
params.BABEL_yaml_dir =  "/shared-workspace/jfisher2/conda/nextflow_envs/BABEL_exported_env-fc3dbdd84f4893d1fca0273cc9f14903"

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
  publishDir "${launchDir}/output/BABEL", mode: 'copy'  // Add this
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

