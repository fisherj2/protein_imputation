
#source relevant python functions


modified_cTPnet=function(data,data_type='Seurat2',model_file_path,dprotein=24, protein_names='', rna_names='',rna_dim=0 ){
  cat('Start data preprocessing...\n')
	if(data_type=='Seurat2'){
		X=preprocess_seurat2(data, dprotein, rna_names)
	}else if(data_type=='Seurat3'){
		X=preprocess_seurat3(data, dprotein, rna_names)
	}else if (data_type=='matrix'|data_type=='dataframe'){
		X=preprocess_matrix(data, dprotein, rna_names)
	}else{
		stop('Error: unrecognizable data_type argument. Need to be one of the four\n
		      options: Seurat2, Seurat3, matrix, dataframe. You can check your \n
		      Seurat version by sessionInfo()\n')
	}
  cat('Start imputation. Running python ...\n')
  #changed this section to run ou modified functions
	#ctpnet <- reticulate::import("ctpnet", convert = F)
	y_pred=modified_predict(X,model_file_path,dprotein,protein_names, rna_dim)
	cat('Postprocess...\n')
	if(data_type=='Seurat2'){
	  data=postprocess_seurat2(data,y_pred)
	}else if(data_type=='Seurat3'){
	  data=postprocess_seurat3(data,y_pred)
	}else{
	  data=postprocess_matrix(y_pred)
	}
	cat('Done!\n')
	return(data)
}



