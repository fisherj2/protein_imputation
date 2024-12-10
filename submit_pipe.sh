#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --job-name=prot_prediction
#SBATCH -e parent_%A_%a.stderr
#SBATCH -o parent_%A_%a.stdout
#SBATCH --partition=Posit-Slurm-m5-32g-8

# ---- setting input_file ----
# running on /scratch/jfisher2/protein_prediction/source_data/Hao_et_al_Seurat_obj.Rdata will use the full, large data
# optionally can change to /scratch/jfisher2/protein_prediction/source_data/Hao_et_al_Seurat_obj_minimal_testing.Rdata  to run on smaller test data

# ----setting profile ----
# profile = standard will run each analysis step one, on one partition into training and test sets. See config file for memory specifications
# profile = benchmarking will run pipeline many times on many different partitions of the data, trying to evaluate performance. See config file for memory specifications
nextflow run /scratch/jfisher2/protein_prediction/main.nf --input_file /scratch/jfisher2/protein_prediction/source_data/Hao_et_al_Seurat_obj.Rdata -profile benchmarking 