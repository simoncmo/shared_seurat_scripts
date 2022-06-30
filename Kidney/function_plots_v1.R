## 
# Function
# 7/1/2021
library('Seurat')
message(str_glue('Function in this file inlcudes: 
	line_plot_list,
	line_plot'))

line_plot_list = function(obj, features, ...){
	map(features, function(feature){
		line_plot(obj, feature, ...)
	}) %>% setNames(features)
}

line_plot =  function(obj, feature, group_by, idents, split_by, split_idents, assay='SCT', x_axis = 'Age'){
	# Add expression to meta
	gene_exp  = GetAssayData(obj, assay= assay, slot = 'data')[feature,]
	exp_meta_df = cbind(obj@meta.data, gene_exp)
	# idents
	if(missing(group_by)) stop('Need group_by argument to plot')
	if(missing(idents)){
		idents = unique(obj@meta.data[[group_by]])
		message(str_glue('Missing idents. Using all {length(idents)} idents'))
	}
	# Filter meta 
	exp_meta_df = exp_meta_df %>% filter(!!as.symbol(all_of(group_by)) %in% idents)
	if(!missing(split_idents)){
		exp_meta_df = exp_meta_df %>% filter(!!as.symbol(all_of(split_by)) %in% split_idents)
	}
	# plot
	p = ggplot(exp_meta_df, aes_string(x = x_axis, y = 'gene_exp')) + 
	#geom_boxplot()+
	#geom_violin()+
	geom_jitter(size = 0.1, aes_string(color = x_axis))+
	stat_summary(fun = 'mean', geom = 'line', aes(group=1), color = 'gray40',, size = 0.3) + 
	stat_summary(fun = 'mean', geom='point', color = 'gray40', size = 2) + 
	stat_summary(fun.data = mean_se, geom = "errorbar", color = 'gray40',, size = 0.3) + 
	labs(title = feature) + 
	ggthemes::theme_clean() + 
	theme(panel.background = element_rect(),
		  panel.grid = element_blank()
			)

	if(!missing(split_by)){
		p = p + facet_grid(cols = vars(!!as.symbol(split_by)),
							rows = vars(!!as.symbol(group_by)),
							scales = 'free_y')
	}else{
		p = p + facet_grid(rows = vars(!!as.symbol(group_by)),
							scales = 'free_y')
	}

	#Plot
	p
}