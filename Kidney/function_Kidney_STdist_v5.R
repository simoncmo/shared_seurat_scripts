library(Seurat)
library(qs)
library(tidyverse)
library(patchwork)

source('/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/Shared_resource/script/function_spatial_plot_helper.R')
make_NMK_plot = function(st_merge, st_male, st_female, gene_to_plot, male_sm_set, male_sm_set_ptsize, female_sm_set, female_sm_set_ptsize){
	########### NMK ########################
	########### individual strips ########################
	# Spatial Dis plot
	p_nmk_M = SpatialFeaturePlotPlus(
		obj = st_male, 
		expression_ref_obj = st_merge,
		features = gene_to_plot, 
		images = male_sm_set, 
		pt.size.factors = male_sm_set_ptsize,
		sample_labels   = c("E165", "P0", "1W", "2W", "3W", "12W", "52W", "113W"),
		feature_position = 'row',
		color_palette = 'viridis',
		assemble=F)

	p_nmk_F = SpatialFeaturePlotPlus(
		obj = st_female, 
		expression_ref_obj = st_merge, 
		features = gene_to_plot, 
		images = female_sm_set, 
		pt.size.factors = female_sm_set_ptsize,
		# sample_labels   = female_sm_set_labels,
		feature_position = 'row',
		color_palette = 'viridis',
		assemble=F)

	p_filler = SpatialFeaturePlotPlus(st_merge,
		features = gene_to_plot,
		images = male_sm_set[[1]],
		pt.size.factors = 0,
		feature_position = 'row',
		color_palette = 'viridis',
		assemble=F
		)

	########### PLOT ASSEMBLE ######################
	############### NMK ############
	p_nmk_MF_gene = lapply(1:length(gene_to_plot), function(i){
		# Male
		M_text          = wrap_elements(grid::textGrob('M'))
		M_layout        = 'TPPPPPPPPPPPPPP'
		M_row_formatted = wrap_plots(p_nmk_M[[i]],nrow=1)
		M_row_formatted = wrap_plots(T=M_text, P=M_row_formatted, design=M_layout, heights = unit(1.2,'inches'))  # Add M on left
		# Female
		F_text          = wrap_elements(grid::textGrob('F'))
		F_layout        = 'TGPPPPPPPPPPPPPP'
		F_filler        = p_filler[[i]][[1]] + theme(legend.position = 'none')
		F_row_formatted = wrap_plots(p_nmk_F[[i]],nrow=1) #& theme(plot.title = element_blank())
		F_row_formatted = wrap_plots(T=F_text,G=F_filler, P=F_row_formatted, design=F_layout, heights = unit(1.2,'inches')) # Add F on left
		F_row_formatted = F_row_formatted & theme(plot.title = element_blank(), legend.position = 'none' )
		# Plot
		M_row_formatted/F_row_formatted
	})
	p_nmk_wrapped = wrap_plots(p_nmk_MF_gene, nrow=length(gene_to_plot)) #& theme(legend.position = 'left')
	return(p_nmk_wrapped)
}

make_AKI_plot = function(st_merge, st_aki, gene_to_plot, aki_male_set, aki_male_set_ptsize, aki_male_set_labels, aki_female_set, aki_female_set_ptsize){
	########### AKI ######################
	# Spatial Dis plot
	p_aki_M = SpatialFeaturePlotPlus(st_aki, 
		features = gene_to_plot, 
		expression_ref_obj = st_merge, # Use Normal as scale reference
		images = aki_male_set, 
		pt.size.factors = aki_male_set_ptsize,
		sample_labels   = aki_male_set_labels,
		feature_position = 'row',
		color_palette = 'viridis',
		assemble=F)
	p_aki_F = SpatialFeaturePlotPlus(st_aki, 
		features = gene_to_plot, 
		expression_ref_obj = st_merge, # Use Normal as scale reference
		images = aki_female_set, 
		pt.size.factors = aki_female_set_ptsize,
		sample_labels   = aki_male_set_labels,
		feature_position = 'row',
		color_palette = 'viridis',
		assemble=F)

	########### PLOT ASSEMBLE ######################
	############ AKI ############
	p_aki_MF_gene = lapply(1:length(p_aki_M), function(i){
		#Add Gender label
		p_aki_M_text = wrap_elements(grid::textGrob('M'))
		p_aki_M[[i]] = wrap_plots(p_aki_M[[i]])
		p_aki_F_text = wrap_elements(grid::textGrob('F'))
		p_aki_F[[i]] = wrap_plots(p_aki_F[[i]])
		#Layout
		design   = "QMMMMMMM
				    WFFFFFFF"
		wrap_plots(Q = p_aki_M_text,
				   M = p_aki_M[[i]], 
				   W = p_aki_F_text, 
				   F = p_aki_F[[i]], design = design)
	})
	return(p_aki_MF_gene[[1]])
}

plot_NMK_AKI  = function(st_merge, st_male, st_female, st_aki, gene_to_plot, title ="", subtitle="", plot_group){
	print(str_glue("Plotting {gene_to_plot}..."))
	# Check if gene exist 
	if(!gene_to_plot %in% rownames(st_merge)){
		message(str_glue("Warning! {gene_to_plot} not found in the obj. Skipping"))
		return()
	}
	##### Preset parameters ############
	today_date = format(Sys.Date(), format = "%m%d%Y")
	##### NMK ##########################
	# Names
	male_sm_set   		  = c("NMK0205E165MU1","NMK0108P0FMU1","NMK0115W1FMU1","NMK0111W2MU1","NMK1201W3MU1","NMK0129W12MU1","NMK1207W52MU1","NMK0226W113MU1")
	female_sm_set 		  = c("NMK0318E165FU1","NMK0108P0FMU1","NMK0115W1FMU1","NMK0111W2FU1","NMK1201W3FU1",'NMK0129W12FU1',"NMK1207W52FU1", 'NMK0430W113FU1')
	# Size
	male_sm_set_ptsize    = c(3, 6.2, 3.2, 1.7 ,1.7 ,2.1 ,2 ,1.5 )
	male_sm_set_labels    = str_c(c("E165", "0W", "1W", "2W", "3W", "12W", "52W", "113W"),"-M")
	female_sm_set_ptsize  = c(3, 4.8, 3.2, 1.7 ,1.7 ,1.8 ,1.6, 2.3)
	female_sm_set_labels  = str_c(c("E165", "0W", "1W", "2W", "3W", "12W", "52W", "113W"),"-F")

	##### AKI ##########################
	# Names
	aki_male_set   		  = c("NMK0129M3MU1", "NMK0129AKID3MU1", "NMK0203AKID8MU1")
	aki_female_set 		  = c("NMK0129M3FU1", "NMK0129AKID3FU1", "NMK1225AKID8FU1")
	# Size
	aki_male_set_ptsize     = c(2, 1.6 ,1.8)
	aki_female_set_ptsize   = c(1.8, 1.8 ,1.8)
	aki_male_set_labels     = c("Control", "3-day", '8-day')

	##### Run and get figures
	# Check before plot
	p_nmk_wrapped = if(gene_to_plot %in% rownames(st_merge)){
		make_NMK_plot(st_merge, st_male, st_female, gene_to_plot, male_sm_set, male_sm_set_ptsize, female_sm_set, female_sm_set_ptsize)
	}else{
		message(str_glue('{gene_to_plot} not found in NMK samples'))
		NULL
	}
	
	p_aki_MF_gene = if(gene_to_plot %in% rownames(st_aki)){
		make_AKI_plot(st_merge, st_aki, gene_to_plot, aki_male_set, aki_male_set_ptsize, aki_male_set_labels, aki_female_set, aki_female_set_ptsize)
	}else{
		message(str_glue('{gene_to_plot} not found in AKI sample'))
		NULL
	}
	
	# Choose plot mode
	if(!is.null(p_nmk_wrapped) & !is.null(p_aki_MF_gene)){
		design = "NNNNNNNNNAAAA"
		p =wrap_plots(N = p_nmk_wrapped, A =p_aki_MF_gene, design = design) + plot_annotation(
				title = title,
				subtitle = subtitle
				)
	}else if(!is.null(p_nmk_wrapped) & is.null(p_aki_MF_gene)){
		message(str_glue("{gene_to_plot} only found in NMK. NMK only mode"))
		design = "NNNNNNNNN"
		p =wrap_plots(N = p_nmk_wrapped, design = design) + plot_annotation(
				title = title,
				subtitle = subtitle
				)
	}else if(is.null(p_nmk_wrapped) & !is.null(p_aki_MF_gene)){
		design = "AAAA"
		p =wrap_plots( A =p_aki_MF_gene, design = design) + plot_annotation(
				title = title,
				subtitle = subtitle
				)
	}else{
		message(str_glue("{gene_to_plot} no found in any object. Skip"))
		return()
	}

	
	##### Set Plot file
	if(!missing(plot_group)){
		plot_folder = str_glue("figure/STExpression/{plot_group}/{today_date}")
		dir.create(plot_folder, recursive=T)
		plot_name = str_glue("{plot_folder}/NMK_normal_AKI_{plot_group}_{gene_to_plot}_{today_date}.pdf")
	}else{
		plot_folder = str_glue("figure/STExpression/Untitle_STplot/{today_date}")
		dir.create(plot_folder, recursive=T)
		plot_name = str_glue("{plot_folder}/NMK_normal_AKI_{gene_to_plot}_{today_date}.pdf")
	}
	
	##### Plot
	pdf(plot_name,width =20, height = 4)
	print(p)
	dev.off()
}


make_expression_plot = function(st_merge, gene_of_interest, split_by, split_filter, error_bar=F){
	# Meta
	if(missing(split_by)){
		line_df_nmk = st_merge@meta.data[,c('Gender','Age_new')]
	}else{
		line_df_nmk = st_merge@meta.data[,c('Gender','Age_new', all_of(split_by))]
	}
	# Expression
	line_df_nmk$Expression = st_merge %>% GetAssayData %>% .[gene_of_interest,] 
	line_df_nmk = filter(line_df_nmk, !is.na(Gender))
	line_df_nmk = mutate(line_df_nmk, Age_new = factor(Age_new, levels = c('E165','P0','W1','W2','W3','W14','W52','W113')))
	# Filter splits
	if(!missing(split_by)&!missing(split_filter)){
		line_df_nmk = filter(line_df_nmk, !!as.symbol(all_of(split_by)) %in% split_filter)
		# line_df_nmk = mutate(line_df_nmk, !!as.symbol(all_of(split_by)) = factor( !!as.symbol(all_of(split_by)), levels = split_filter  ) )
	}

	# Summarize
	if(missing(split_by)){
		line_df_nmk_sum = line_df_nmk %>% group_by(Gender, Age_new) %>% summarize(exp_mean = mean(Expression), exp_sd = sd(Expression))
	}else{
		line_df_nmk_sum = line_df_nmk %>% group_by(Gender, Age_new, !!as.symbol(all_of(split_by))) %>% summarize(exp_mean = mean(Expression), exp_sd = sd(Expression))
	}
	
	# plot
	p = ggplot(line_df_nmk_sum, aes(x = Age_new, y= exp_mean, color = Gender, group=Gender)) + 
			geom_line()+
			geom_point()+
			scale_color_manual(values = c('M'='dodgerblue','F'='brown1'),)+
			theme_minimal()+
			labs(x = 'Age', y = gene_of_interest)+
			theme(panel.grid = element_blank(),
				axis.line = element_line())
			
	if(error_bar){
		p = p + geom_errorbar(aes(ymin=exp_mean-exp_sd, ymax=exp_mean+exp_sd), width=.2,
                 position=position_dodge(0.05)) 
	}

	if(!missing(split_by)){
		p = p + facet_wrap(reformulate(split_by,"."), scales = 'free')
	}
	p
}
