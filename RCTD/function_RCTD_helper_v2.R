## RCTD helper function
## 06/03/2021
library(tidyverse)
library(RCTD)


RCTD_obj2list = function(RCTD_result_obj){
	message("This method drop the reference object to reduce size")
	RCTD_result_list = list()
	RCTD_result_list$spatialRNA     = RCTD_result_obj@spatialRNA
	RCTD_result_list$config         = RCTD_result_obj@config
	RCTD_result_list$cell_type_info = RCTD_result_obj@cell_type_info
	RCTD_result_list$internal_vars  = RCTD_result_obj@internal_vars
	RCTD_result_list$results        = RCTD_result_obj@results
	# normalize the cell type proportions to sum to 1.
	RCTD_result_list$norm_weights = sweep(RCTD_result_obj@results$weights, 1, rowSums(RCTD_result_obj@results$weights), '/') 
	message(str_glue("Original size: {format(object.size(RCTD_result_obj),units='Mb')}"))
	message(str_glue("Reduced size: {format(object.size(RCTD_result_list),units='Mb')}"))
	structure(RCTD_result_list, class = 'RCTDshrink')
}

getRCTD_result = function(st_obj, snrna_obj, colu_cell_type='cell_type_Abbr', mode ='full',cell_count_filter = 25){
	# filter out cell type <= 25
	message(str_glue('filtering out cell type <= {cell_count_filter}'))
	celltype_to_keep  = names(table(snrna_obj[[colu_cell_type]]))[table(snrna_obj[[colu_cell_type]])>cell_count_filter]
	Idents(snrna_obj) = colu_cell_type
	snrna_obj         = subset(snrna_obj, idents = celltype_to_keep)
    # Clean factor levels
    snrna_obj@meta.data[[colu_cell_type]] = snrna_obj@meta.data[[colu_cell_type]] %>% 
                                            factor(., levels = unique(.))
	# Create ref obj
	sn_counts     = GetAssayData(object = snrna_obj, slot = "counts") # Get Raw Count
	sn_meta       = snrna_obj@meta.data
	sn_cell_types = sn_meta[[colu_cell_type]]   %>% setNames(rownames(sn_meta)) %>% as.factor(.)
	sn_nUMI       = colSums(sn_counts)           # This is optional. create nUMI named list
	snrna_ref <- Reference(sn_counts, sn_cell_types, sn_nUMI)

	## ST construnction
	st_counts = GetAssayData(st_obj, slot = 'counts')
	st_coords = map(unique(st_obj$orig.ident), ~ GetTissueCoordinates(st_obj, image =.x)) %>% 
			bind_rows() %>%
			setNames(c("x","y")) # Get All the image rows in the obj
	st_nUMI   = colSums(st_counts) # Optional. In this case, total counts per pixel is nUMI

	### Create SpatialRNA object
	st_puck <- SpatialRNA(st_coords, st_counts, st_nUMI)

	### Create RCTD obj
	print(Sys.time()-ptm);ptm =Sys.time()
	RCTD_obj <- create.RCTD(st_puck, snrna_ref, max_cores = 1)

	### RUN RCTD
	print(Sys.time()-ptm);ptm =Sys.time(); print('Starting RCTD')
	RCTD_obj <- run.RCTD(RCTD_obj, doublet_mode = mode)
	print(Sys.time()-ptm);ptm =Sys.time(); print('Done!')

	RCTD_list = RCTD_obj2list(RCTD_obj)
	return(RCTD_list)
}

### RCTD result to Assay and Meta 
## 10/18/2021
RCTD2assay = function(x, meta, ...){
	UseMethod('RCTD2assay')
}

RCTD2assay.RCTD = function(x, meta, ...){
    x = x@results$weights %>% as.matrix
    NextMethod('RCTD2assay')
}
RCTD2assay.dgCMatrix = function(x, meta, ...){
    x = as.matrix(x) 
    NextMethod('RCTD2assay')
}

RCTD2assay.matrix = function(x, meta, ...){
    x = as.data.frame(x) 
    NextMethod('RCTD2assay')
}


RCTD2assay.RCTDshrink = function(x, meta, ...){
	message('1 shrinked RCTD results')
	x = x$norm_weights %>% as.data.frame
	NextMethod('RCTD2assay')
}
RCTD2assay.list = function(x, add_id = T, ...){
	message('List of RCTD results')
    x = imap(x, function(result, name){
        df = result$norm_weights %>% as.data.frame
        if(add_id){
            message("Add item name to id. Change this by setting add_id = F")
            rownames(df) = str_c(name, rownames(df), sep='_')
        }
        df
    }) %>% bind_rows()
	NextMethod('RCTD2assay')
}
RCTD2assay.default = function(x, meta, mode = c('all', 'intersect'), na_fill = 0, return_as = c('matrix','assay'), sample_name,...){
	# Merge mode
	mode    = match.arg(mode)
	rowFn   = switch(mode, 'intersect' = intersect, 'all' = union)
    
    # rowname prefix
    if(!missing(sample_name)) rownames(x) = str_c(sample_name, rownames(x), sep='_')

	# Merge
	if(!missing(meta)){
		rowName = rowFn(rownames(meta), rownames(x))
		x = x[rowName, ] 
		rownames(x) = rowName
	}
	# Fill na
	x[is.na(x)] = na_fill

	# Return
	if(match.arg(return_as) =='assay'){ return(CreateAssayObject(t(x)))
	}else{ return(x) }


}

### Meta
ratio2meta = function(x, col_name= 'RCTD'){
	df = data.frame(row.names = rownames(x),
					name   = colnames(x)[apply(x, 1, function(row) ifelse(length(unique(row))==1, NA, which.max(row) ))],  # All NA or all 0
					weight = apply(x, 1, function(row) ifelse(all(is.na(row)), NA, max(row, na.rm=T)))
					)
	df = setNames(df, paste(col_name, names(df), sep='_'))
	df
}

RCTD2meta = function(x, ...){
	UseMethod('RCTD2meta')
}
RCTD2meta.list = function(x, mode = 'all', ...){
	x = RCTD2assay(x, mode = 'all', ...)
	NextMethod('RCTD2meta')
}
RCTD2meta.RCTDshrink = function(x, ...){
	x = RCTD2assay(x, mode = 'all', ...)
	NextMethod('RCTD2meta')
}
RCTD2meta.Assay = function(x,...){
	x = GetAssayData(x) %>% as.matrix %>% t
	NextMethod('RCTD2meta')
}
RCTD2meta.default = function(x, col_name, ...){
	ratio2meta(x, col_name, ...)
}

                                               
# 12/13/2021
fix_cellname = function(names){ str_replace_all(names,'-','_')}

fix_cellname_mtx = function(obj, mode = c('column','row')){
    mode = match.arg(mode)
    if(mode == 'column'){
        colnames(obj) = colnames(obj) %>% fix_cellname()
    }else{
        rownames(obj) = rownames(obj) %>% fix_cellname()
    }
    obj
}

integrate_RCTDassay_mtx = function(orig_mtx, new_mtx, fix_name = T){
if('Assay' %in% class(orig_mtx)) orig_mtx = GetAssayData(orig_mtx) %>% as.matrix %>% t 
if('Assay' %in% class(new_mtx))  new_mtx = GetAssayData(new_mtx) %>% as.matrix %>% t 
if(fix_name){
    print("All column id name FIXED. Double check result if this is desired")
    orig_mtx = orig_mtx %>% fix_cellname_mtx(mode = 'column')
    new_mtx  = new_mtx %>% fix_cellname_mtx(mode = 'column')
}
    all_ids  = union(rownames(orig_mtx), rownames(new_mtx))
    all_cols = union(colnames(orig_mtx), colnames(new_mtx))
    all_mtx  = matrix(ncol = length(all_cols), nrow = length(all_ids)) %>% as.data.frame
    colnames(all_mtx) = all_cols
    rownames(all_mtx) = all_ids
    all_mtx[rownames(orig_mtx), colnames(orig_mtx)] = orig_mtx
    all_mtx[rownames(new_mtx), colnames(new_mtx)]   = new_mtx
    all_mtx[is.na(all_mtx)] = 0 #fill na 0
    return(all_mtx)
}

integrate_RCTDmeta = function(orig_meta, new_meta, target_cols, create_new_col=T){
    if(create_new_col){
        new_cols   = setdiff(colnames(new_meta), colnames(orig_meta))
        new_col_df = data.frame(matrix(data=NA, ncol =length(new_cols), nrow=nrow(orig_meta)), 
                                row.names = rownames(orig_meta)) %>% 
                    setNames(new_cols) # NEW NA df, same rownames as orig_meta
        orig_meta = cbind(orig_meta, new_col_df[rownames(orig_meta), ]) # ADD to orig_meta
                #print(head(orig_meta))
    }
    if(missing(target_cols)){
        target_cols = intersect(colnames(orig_meta), colnames(new_meta))
        print(str_glue("Add new values to {toString(target_cols)}"))
    }
    new_ids  = rownames(new_meta)
    orig_meta[new_ids, target_cols] = new_meta[new_ids, target_cols]
    return(orig_meta)
}

order_mtx_ids = function(mtx, ordered_ids){
    if(!all(ordered_ids %in% rownames(mtx))) stop("not all ordered_ids in the mtx. Please check")
    mtx[ordered_ids, ]
}
                                               
integrate_RCTDmeta_to_obj = function(obj, new_meta){
	new_whole_meta = integrate_RCTDmeta(obj@meta.data, new_meta)
	obj@meta.data = new_whole_meta[rownames(obj@meta.data), ]
	obj
}