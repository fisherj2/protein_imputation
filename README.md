# protein_imputation
Repository holding the code associated with the work 'Machine Learning Predictions Surpass Individual mRNAs as a Proxy of Single-cell Protein Expression'

We tested the performance of several protein imputation methods on a selection of multimodal datasets:

* [Hao et al. PBMC CITEseq data](https://atlas.fredhutch.org/nygc/multimodal-pbmc)
* [Stephenson et al. PBMC CITEseq data](https://pubmed.ncbi.nlm.nih.gov/33879890/)
* [Leader et al. NSCLC CITEseq data](https://pubmed.ncbi.nlm.nih.gov/34767762/)
* [Zhang et al. RA synvoium CITEseq data](https://www.nature.com/articles/s41590-024-01782-4)

Our analyses were implemented using a Nextflow pipeline. You can find more details regarding Nextflow development and execution in the [documentation](https://www.nextflow.io/docs/latest/index.html). To understand how we executed the pipeline, consult the submit_pipe.sh file. 

Each protein imputation method was implemented using a separate conda environment, to avoid version conflicts. These environments are defined in yml files alongside the method execution scripts in the modules folder. 