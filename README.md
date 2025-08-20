# protein_imputation
Repository holding the code associated with the work 'Machine Learning Predictions Surpass Individual mRNAs as a Proxy of Single-cell Protein Expression'

We tested the performance of several protein imputation methods on the [Hao et al PBMC CITEseq data](https://atlas.fredhutch.org/nygc/multimodal-pbmc). For ease of repetition our analysis was implemented using a Nextflow pipeline. You can find more details regarding Nextflow development and execution in the [documentation](https://www.nextflow.io/docs/latest/index.html). To understand how we executed the pipeline, consult the submit_pipe.sh file. 

Each protein imputation method was implemented using a separate conda environment, to avoid version conflicts. These environments are defined in yml files alongside the method execution scripts in the modules folder. 