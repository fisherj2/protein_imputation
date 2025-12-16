
params.scTranslator_yaml_dir =  "${moduleDir}/env/scTranslator_env_exported.yml"

// NOTE: encountered issues when trusting nextflow to build the environment itself
// need to build beforehand, test, and then link it here. see scTranslator docs for guidance.
params.scTranslator_env_dir = "/shared-workspace/jfisher2/conda/envs/performer"
process scTranslator_init_conda{
  label 'small_mem'
  if ( params.do_scTranslator) {
    //path to desired conda environment
    conda params.scTranslator_yaml_dir 
  }

  input:
   stdin

  output:
  stdout


  script: 
  """

  source ${moduleDir}/bin/check_cuda_and_reinstall.sh
  
  """
  
}

//this process defines scipenn method training
process scTranslator_build_h5ad{
  label 'medium_mem'
  //path to desired conda environment
  conda params.scTranslator_env_dir
  publishDir "${launchDir}/output/scTranslator", mode: 'copy'  
  
  
  input:

      val file 

  output:

      path "*.h5ad", emit: h5ad_files
  


  script: 
  """
    echo "building h5ad objects"
    eval "python  ${moduleDir}/bin/build_h5ad.py --basedir $projectDir  --bench ${params.dobenchmark} --files  '${file}' --moduledir=${moduleDir}"
  
  """
  
}

process scTranslator_map_h5ad{
  label 'medium_mem'
  //path to desired conda environment
  conda params.scTranslator_env_dir
  publishDir "${launchDir}/output/scTranslator", mode: 'copy'  
  
  input:

      val built_file

  output:

      path "*_mapped.h5ad"
  


  script: 
  """
python ${moduleDir}/bin/code/model/data_preprocessing_ID_convert.py \
--origin_gene_type='human_gene_symbol' \
--origin_gene_column='index' \
--data_path=${built_file} \
--moduledir='${moduleDir}'
  """
  
}


//this process defines scipenn method training
process scTranslator_train{
  label 'gpu'
  //path to desired conda environment
  conda params.scTranslator_env_dir
  publishDir "${launchDir}/output/scTranslator", mode: 'copy'  

  input:

      val file 

  output:

      path "*.csv"
      path "*.pt"

  script: 
  """
  #isolate the h5ad files
  
  #pass them to conversion and save to path
  
  #pass those converted files to the finetuning function

  eval "python -m torch.distributed.launch --nnodes=1 --node_rank=0 --nproc_per_node 1 --master_port 23333 \
${moduleDir}/bin/predictADT_scTranslator.py --epoch=100 --frac_finetune_test=0.1 --fix_set  --batch_size=4 \
--pretrain_checkpoint='/home/jfisher2/workspace/analyses/protein_prediction_publication/repo/modules/local/scTranslator/checkpoint/scTranslator_2M.pt' \
--RNA_path='' \
--Pro_path='' \
--basedir $projectDir \
--moduledir $moduleDir \
--bench ${params.dobenchmark} \
--files  '${file}'
  "


  """
}

process scTranslator_fewshot{
  label 'gpu'
  //path to desired conda environment
  conda params.scTranslator_env_dir
  publishDir "${launchDir}/output/scTranslator", mode: 'copy'  

  input:

      val file 

  output:

      path "*.csv"
      path "*.pt"

  script: 
  """

  eval "python -m torch.distributed.launch --nnodes=1 --node_rank=0 --nproc_per_node 1 --master_port 23333 \
${moduleDir}/bin/predictADT_scTranslator_fewshot.py --epoch=100 --frac_finetune_test=0.1 --fix_set  --batch_size=4 \
--pretrain_checkpoint='/home/jfisher2/workspace/analyses/protein_prediction_publication/repo/modules/local/scTranslator/checkpoint/scTranslator_2M.pt' \
--RNA_path='' \
--Pro_path='' \
--basedir $projectDir \
--moduledir $moduleDir \
--bench ${params.dobenchmark} \
--files  '${file}' \
--ncells 1000
  "


  """
}


process scTranslator_notrain{
  label 'gpu'
  //path to desired conda environment
  conda params.scTranslator_env_dir
  publishDir "${launchDir}/output/scTranslator", mode: 'copy'  

  input:

      val file 

  output:

      path "*.csv"


  script: 
  """

  #pass those converted files to the finetuning function

  eval "python ${moduleDir}/bin/predictADT_scTranslator_no_finetune.py \
--pretrain_checkpoint='/home/jfisher2/workspace/analyses/protein_prediction_publication/repo/modules/local/scTranslator/checkpoint/scTranslator_2M.pt' \
--basedir $projectDir \
--moduledir $moduleDir \
--fix_set \
--bench ${params.dobenchmark} \
--files  '${file}'
  "


  """
}