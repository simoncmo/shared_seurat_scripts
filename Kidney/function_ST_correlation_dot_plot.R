# 08042021

get_all_Agesex_cor_df = function(obj, group_by = 'Age_sex', ...){
	# Add Age sex meta
	if(!'Age_sex' %in% names(obj@meta.data)) obj$Age_sex = str_c(obj$Age, obj$Sex, sep='_')

	# Remove 'P0_NA'
	age_sex_groups = setdiff(unique(obj$Age_sex), 'P0_NA')
	
	# Make List
	corr_df_list = map(age_sex_groups, function(age_sex){
		get_cor_df(st_merge, ident = age_sex, complete_matrix=T, ...)
	}) %>% setNames(age_sex_groups)

	# Remove NA groups
	corr_df_list = corr_df_list[unlist(map(corr_df_list, ~typeof(.) =='list'))]

	# Combine
	corr_df_all = bind_rows(corr_df_list) %>% dplyr::rename(feature1= row, feature1 = col)
	# Fix age order
	mutate(corr_df_all, Age = factor(Age, levels = c('E165','P0','W1','W2','W3','W12','W52','W113')))

}

get_celltype_cor_df = function(obj, ident, group_by = 'Age_sex', features ,assay_name = 'RCTD_doubletRM_v3', test_sample_threshold = 5, complete_matrix =T, cor_method='complete.obs', split_col_by = '_'){
	get_cor_df(obj=obj, ident=ident, group_by=group_by, features=features, assay_name=assay_name,
			   test_sample_threshold=test_sample_threshold,
			   complete_matrix=complete_matrix,
			   cor_method=cor_method,
			   split_col_by=split_col_by
			)
}
get_exp_cor_df = function(obj, ident, group_by = 'Age_sex', features ,assay_name = 'SCT', test_sample_threshold = 5, complete_matrix =T, cor_method='complete.obs', split_col_by = '_'){
	get_cor_df(obj=obj, ident=ident, group_by=group_by, features=features, assay_name=assay_name,
			   test_sample_threshold=test_sample_threshold,
			   complete_matrix=complete_matrix,
			   cor_method=cor_method,
			   split_col_by=split_col_by
			)
}
get_cor_df = function(obj, ident, group_by = 'Age_sex', features ,assay_name = 'SCT', test_sample_threshold = 5, complete_matrix =T, cor_method='complete.obs', split_col_by = '_'){
	# Check
	if(missing(ident)) print('Missing ident. Default group_by = Age_sex. Ident example: W113_F')
	# Select cells
	cell_ids = filter(obj@meta.data, !!as.symbol(group_by) %in% ident) %>% rownames
	if(length(cell_ids) < test_sample_threshold){
		message(str_glue("Test group {ident} has {length(cell_ids)} samples, less than {test_sample_threshold}. Skip."))
		return(NA)
	}

	# Select features
	if(missing(features)){
		features = rownames(GetAssay(obj, assay = assay_name))
		if(length(features)>1000){
		stop('Too many (>1000) features selected. Might wanna provide "features" to compare instead?')
	}
	message(str_glue('Missing features to compare, use all {length(features)} features from {assay_name}'))
	}
	
	# Get Data	
	DefaultAssay(obj) = assay_name
	cell_percent_matrix = FetchData(obj, vars = features, cells = cell_ids)

	# Remove NA only cell type/column
	non_na_cell = apply(cell_percent_matrix, 2, function(x) sum(!is.na(x))!=0)
	cell_percent_matrix = cell_percent_matrix[,non_na_cell]

	message(unique(apply(cell_percent_matrix, 2, function(x) sum(!is.na(x)))))

	# # Correlation, p dataframe
	cor_p_df = matrix2cor_df(cell_percent_matrix ,complete_matrix =complete_matrix, cor_method=cor_method)

	# Add meta
	cor_p_df[[group_by]] = ident

	# split column
	if(str_detect(ident, split_col_by)){
		into_cols = str_split(group_by, pattern = split_col_by) %>% unlist
		cor_p_df  = cor_p_df %>% separate( (!!sym(group_by)), into= into_cols, sep='_') 
	}
	
	# rename and output
	cor_p_df %>% dplyr::rename(corr = r )

}

triangle_to_df = function(m, value_name ='corr'){
	data.frame(row=rownames(m)[row(m)[upper.tri(m)]], 
           col=colnames(m)[col(m)[upper.tri(m)]], 
           corr=m[upper.tri(m)]) %>% setNames(c('row','col',value_name))
}

matrix2cor_df = function(value_matrix, complete_matrix =T, cor_method='complete.obs'){
	# Correlation
	cor_matrix    = cor(value_matrix, use = cor_method)
	cor_df        = triangle_to_df(cor_matrix, value_name = 'r')
	if(complete_matrix){
		cor_df_bottom = cor_df %>% dplyr::rename(row = col, col = row)
		cor_df = bind_rows(cor_df, cor_df_bottom)
	}

	# p
	p_matrix = corrplot::cor.mtest(mat = value_matrix, conf.level = 0.95)
	rownames(p_matrix$p) = rownames(cor_matrix)
	colnames(p_matrix$p) = colnames(cor_matrix)
	p_df = triangle_to_df(p_matrix$p, value_name = 'p')
	if(complete_matrix){
		p_df_bottom = p_df %>% dplyr::rename(row = col, col = row)
		p_df = bind_rows(p_df, p_df_bottom)
	}
	p_df = p_df %>% mutate(fdr = p.adjust(p, method = 'fdr'))
	
	# combine p, cor
	return(merge(cor_df, p_df, by=c('row','col')))
}
