
params.Renv_yaml_dir=  "${moduleDir}/env/default_Renv.yml"

process Renv_init_conda{
  
  label 'small_mem'

  //path to desired conda environment
  conda params.Renv_yaml_dir 


  input:


  output:
  stdout


  script: 
  """
   echo "building default R_env"
  """
  
}


process seurat_anchors{
  label 'big_mem'
  //path to desired conda environment
  conda params.Renv_yaml_dir 
  
  
  input:

      val file 

  output:

      path "*_prediction.csv"

  script: 
  """
   eval "Rscript ${moduleDir}/bin/seurat_anchors.R  $projectDir $launchDir  ${params.dobenchmark} '${file}'"
  """
}
