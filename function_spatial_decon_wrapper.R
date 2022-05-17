# Spatial deconvolution/integration helper
library(Seurat)
library(RCTD)



# Functions version
st_decon_seurat = function(st_merge, snrna, st_transfer_assay, st_transfer_col, sn_cell_type_col){
	##### Label transfer/ integrations ######################################################################
	### Label transfer using Seurat TransferData method -----------
	# Using all 3000 features
	print('Seurat integration transfer Started')
	transfer_anchor   <- FindTransferAnchors(reference = snrna, query = st_merge, normalization.method = "SCT")
	predictions.assay <- TransferData(
	    anchorset = transfer_anchor, 
	    refdata = snrna@meta.data[, sn_cell_type_col], 
	    prediction.assay = TRUE, 
	    weight.reduction = st_merge[["pca"]], 
	    dims = 1:30)
	st_merge[[st_transfer_assay]] = predictions.assay
	st_merge[[st_transfer_col]]   = getDominantAssignment(st_merge, assay_name=st_transfer_assay)
	print('Seurat integration transfer completed.')

	return(st_merge)
}

st_decon_RCTD = function(st_merge, snrna, st_RCTD_assay, st_RCTD_col, sn_cell_type_col, min_cell_type_count=25, st_RCTD_obj, doublet_mode = 'multi'){
	# # --------------------------------------------
	# ### Label transfer with RCTD
	# Prep PCTD ------------
	# Create ref and puck objects
	print('RCTD Deconvolution Started')
	ref_sn  = seurat_to_RCTD_ref_v2(snrna, cell_type_col=sn_cell_type_col, filter_cell_type =T, min_cell_type_count =25)
	puck_st = seurat_to_RCTD_st(st_merge)

	# RUN RCTD - snRNA
	st_RCTD_sn = create.RCTD(puck_st, ref_sn, max_cores = 60, fc_cutoff = 0.3, CELL_MIN_INSTANCE = min_cell_type_count)
	st_RCTD_sn = run.RCTD(st_RCTD_sn, doublet_mode = 'full') ### This step takes hours to run

	# Save RCTD result
	dir.create('obj')
	if(missing(st_RCTD_obj)){st_RCTD_obj = str_glue('RCTD_output_{format(Sys.time(), "%H%M%S")}')}
	qsave(st_RCTD_sn, str_glue('obj/{st_RCTD_obj}.qs'))

	# Save assay to object
	st_merge 					= get_RCTD_proportion(st_RCTD_sn, st_merge, add_to_st=T, assay_name = st_RCTD_assay) #need rerun
	st_merge[[st_RCTD_col]]     = getDominantAssignment(st_merge, assay_name=st_RCTD_assay) 
	print('RCTD Deconvolution completed.')

	return(st_merge)
}

st_decon_prefly_check = function(st_merge, snrna, sn_cell_type_col){
	print('Prefly check ...')
	# Reference column name check
	if(sn_cell_type_col %in% names(snrna@meta.data)){
		print( str_glue("found {length(snrna@meta.data[, sn_cell_type_col])} cells."))
		print( str_glue("found {toString(unique(snrna@meta.data[, sn_cell_type_col]))} of total
			{length(unique(snrna@meta.data[, sn_cell_type_col]))} unique annotations."))
	}else{
		print(str_glue('ERROR: column "{sn_cell_type_col}" not found in the reference object.'))
		print(str_glue('Available columns in the reference are {toString(colnames(snrna@meta.data))}'))
		stop(str_glue('Please choose the right column from the names above.'))
	}

	# Reference size check 
	print(str_glue("Reference obj dimension = {toString(dim(snrna))}"))
	if(nrow(snrna) == 0 | ncol(snrna) == 0){
		stop("Reference size error. Please double check the reference input")
	}

	# ST size check 
	print(str_glue("ST obj dimension = {toString(dim(st_merge))}"))
	if(nrow(st_merge) == 0 | ncol(st_merge) == 0){
		stop("ST size error. Please double check the ST object input")
	}

}

st_decon_saving = function(st_merge, save_type, st_out_obj){
	dir.create('obj/')
	if(save_type == 'qs' & ("package:qs" %in% search())){
		print('Saving st obj as qs file')
		qsave(st_merge, str_glue('obj/{st_out_obj}.qs'))
	}else{
		print('Saving st obj as RDS file')
		saveRDS(st_merge, str_glue('obj/{st_out_obj}.rds'))
	}
}

st_decon_wrapper = function(st_merge, snrna, snrna_assay_name = "snrna", sn_cell_type_col, st_run_version = "v1", save_type = 'rds',
	decon_methods = c('seurat','RCTD'), min_cell_type_count=25, RCTD_doublet_mode = 'full'){
	ptm = Sys.time()
	# Set names 
	st_transfer_assay = str_c(snrna_assay_name, 'integration_assay', sep='_')
	st_transfer_col   = str_c(snrna_assay_name, 'integration_max', sep='_')
	st_RCTD_assay	  = str_c(snrna_assay_name, 'RCTD_assay', sep='_')
	st_RCTD_col 	  = str_c(snrna_assay_name, 'RCTD_max', sep='_')
	st_RCTD_obj       = str_c(snrna_assay_name, 'RCTD_obj', sep='_')
	st_out_obj        = str_c('NMK_normal',snrna_assay_name, st_run_version, sep='_')

	# Sanity check
	st_decon_prefly_check(st_merge, snrna, sn_cell_type_col)
	print(Sys.time() - ptm); ptm = Sys.time()

	# Method check
	available_methods = c('seurat','RCTD')
	if(length(intersect(decon_methods,available_methods))==0){
		stop(str_glue("Currenly available methods are {toString decon_methods} only"))
	}

	#Run decon method
	# 1. Seurat integration
	print(Sys.time() - ptm); ptm = Sys.time()
	if('seurat' %in% decon_methods){
		st_merge = st_decon_seurat(st_merge, snrna, st_transfer_assay, st_transfer_col, sn_cell_type_col)
		print(Sys.time() - ptm); ptm = Sys.time()
	}
	if('RCTD' %in% decon_methods){
		# 2. RCTD decon
		# (st_merge, snrna, st_RCTD_assay, st_RCTD_col, sn_cell_type_col, min_cell_type_count=25, st_RCTD_obj)
		st_merge = st_decon_RCTD(st_merge, snrna, st_RCTD_assay, st_RCTD_col, sn_cell_type_col, min_cell_type_count = min_cell_type_count, st_RCTD_obj = st_RCTD_obj, doublet_mode=RCTD_doublet_mode)
		print(Sys.time() - ptm); ptm = Sys.time()
	}
	else{
		stop(str_glue("Currenly available methods are {toString decon_methods} only"))
	}
	
	

	print('All decon completed. Saving object..')
	st_decon_saving(st_merge, save_type, st_out_obj)
	print(Sys.time() - ptm); ptm = Sys.time()
	print('Save completed.')
	return(st_merge)

}

st_decon_wrapper_paired = function(st_merge, snrna, snrna_assay_name = "snrna", sn_cell_type_col, st_run_version = "v1", save_type = 'rds',
	decon_methods = c('seurat','RCTD'), st_split_column_name, snrna_split_column_name){
	#Check columns
	if(missing(st_split_column_name) || missing(snrna_split_column_name)){
		stop("Need both st_split_column_name, and snrna_split_column_name to run!")
	}

	# Check pairs
	paired_groups = intersect(unique(st_merge@meta.data[[st_split_column_name]]), 
								unique(snrna@meta.data[[snrna_split_column_name]]))
	if(length(paired_groups)==0){
		print(str_glue("ST groups = {toString(unique(st_merge@meta.data[[st_split_column_name]])}"))
		print(str_glue("reference groups = {toString(unique(snrna@meta.data[[snrna_split_column_name]])}"))
		stop("No paired groups found in ST and reference. Make sure used the right column.")
	}

	#### Split objects
	## Age groups: 6 groups Embryo(E165), Newborn(P0), Young(1W,2W,3W), Adult(W12-14, M3), Mid Age(W52), Old(W113, M26)
	## snrna
	print("Splitting snRNA object")
	snrna_list <- SplitObject(snrna, split.by = snrna_split_column_name)

	## ST
	print("Splitting ST object")
	st_list <- SplitObject(st_merge, split.by = st_split_column_name)

	# Match deconvolution runs
	st_decon_list = list()
	for(pair in intersect(names(st_list), names(snrna_list))){
		print(str_glue('Running {pair} now ..'))
		st_decon_list[[pair]] = st_decon_wrapper(
		st_merge = st_list[[pair]], 
		snrna    = snrna_list[[pair]], 
		snrna_assay_name = str_glue("{snrna_assay_name}_{pair}"),
		sn_cell_type_col = sn_cell_type_col, 
		st_run_version 	= st_run_version, 
		save_type 		= save_type,
		decon_methods 	= decon_methods)
	}
	print("All pair calculation completed!")
	return(st_decon_list)
}