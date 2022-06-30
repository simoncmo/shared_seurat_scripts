# Handle missing median_umi
fix_median_umi = function(SCTModel_obj){
	err_message = ''
	tryCatch({ test <- validObject(SCTModel_obj) }, 
		error = function(error_message) {
			err_message <<- as.character(error_message)
	})
	missing_medium_umi = str_detect(err_message, 'median_umi')

	if(missing_medium_umi){
		message('Missing medium_umi, calculate again from cell.attributes$umi')
		slot(SCTModel_obj, 'median_umi') = median(SCTModel_obj@cell.attributes$umi)
	}
	return(SCTModel_obj)
}

# Cleaning empty objects
# General purpose
clean_seurat_obj_list = function(obj_list, attirbute_to_check){
	if(missing(attirbute_to_check)) {stop("Need attributes to check for cleaning")}
	# Object type
	obj_type = class(obj_list[[1]])[[1]]

	# Count
	obj_size = map(obj_list, function(object){
		nrow(slot(object, attirbute_to_check)) 
	}) %>% unlist 

	# Remove empty
	if(length(obj_size ==0)  != 0 ){
		message(str_glue('Removing {length(obj_size ==0)} empty object from the {obj_type} object list'))
		obj_list = obj_list[obj_size!=0]
		message(str_glue('{length(obj_list)} {obj_type} object(s) left'))
	}
	obj_list
}

# for SCTModel.list slot
clean_seurat_SCTModel_list = function(sct_model_list){
	clean_seurat_obj_list(obj_list = sct_model_list, attirbute_to_check = 'cell.attributes')
}

# For images slot
clean_seurat_image_list = function(img_list){
	clean_seurat_obj_list(obj_list = img_list, attirbute_to_check = 'coordinates')
}

###########################
### API: 
###########################
# API to clean empty SCTModel and add missing median
FixSeuratSCT = function(obj){
    fix_seurat_SCT(obj)
}
FixSeuratImage = function(obj){
    fix_seurat_image(obj)
}

fix_seurat_SCT = function(obj){
	# Check first 
	if(!'SCT' %in% Assays(obj)){
		message('SCT assay not found. Nothing to fix')
		return(obj)
	}

	# Model list
	sct_model_list = obj$SCT@SCTModel.list
	# 1. clean SCTModel list
	sct_model_list = clean_seurat_SCTModel_list(sct_model_list)

	# 2. fix missing median_umi
	sct_model_list = map(sct_model_list, function(sct_model){
		fix_median_umi(sct_model)	
	})
	
	# Add back and retrun
	obj$SCT@SCTModel.list = sct_model_list

	return(obj)
}

# API to clean empty Suerat image
fix_seurat_image = function(obj){
	# Check first 
	if(is.null(Images(obj))){
		message('No Images in this obj. Not a ST object')
		return(obj)
	}

	# clean image
	obj@images = clean_seurat_image_list(obj@images)
	return(obj)
}