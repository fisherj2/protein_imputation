

params.SPECK_yaml_dir = "${moduleDir}/env/SPECK_env.yml"


process SPECK_init_conda{
  label 'small_mem'
  if ( params.do_speck) {
    //path to desired conda environment
    conda params.SPECK_yaml_dir
  }


  input:
   stdin

  output:
    stdout


  script: 
  """
   echo "building SPECK_env"

  """
  
}



//this process defines SPECK prediction based only on RNA
process SPECK{
  label 'medium_mem'
  //path to desired conda environment
  conda params.SPECK_yaml_dir 
  
  input:

      val file 

  output:

      path "*_prediction.csv"

  script: 
  """
   eval "Rscript ${moduleDir}/bin/SPECK.R  $projectDir ${params.dobenchmark}  '${file}'"
  """
}

