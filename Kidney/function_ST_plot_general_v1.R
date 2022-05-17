## ST genera plotting function ###
## 06/09/2021 ####
## Simon Mo ######
library(Seurat)
#library(scatterpie)
library(RColorBrewer)
####################################
######### Data processing ##########
st_get_expression = function(coord_df, obj, features){
  exp_df = GetAssayData(obj)[features, rownames(coord_df)]
  if(length(features) == 1){
    coord_df[[features]] = exp_df
    coord_df
  }else{ 
    exp_df = exp_df %>% as.matrix %>% t
    cbind(coord_df, exp_df) 
  }
}

st_get_assay_value = function(coord_df, obj, assay_use, features){ # 7/27/2021
	assay_df = GetAssay(obj, assay = assay_use) %>% GetAssayData()
	
	if(!features %in% rownames(assay_df)){
		stop(str_glue('Feature {features} not found in assay {assay_use}'))
	}

	assay_df = assay_df[features, rownames(coord_df)]

  if(length(features) == 1){
    coord_df[[features]] = assay_df
  }else{ 
    assay_df = assay_df %>% as.matrix %>% t
    coord_df = cbind(coord_df, assay_df) 
  }

  # Fix assay name if have special charactor that'd cause issue
  #colnames(coord_df) = make.names(colnames(coord_df))
  return(coord_df)
}

st_add_meta = function(coord_df, obj, features){
  if(intersect(features, colnames(obj@meta.data))==0) stop('Cant find any features in the meta data of obj')
  exp_df = obj@meta.data[rownames(coord_df), features]
  if(length(features) == 1){
    coord_df[[features]] = exp_df
    coord_df
  }else{ 
    exp_df = exp_df %>% as.matrix %>% t
    cbind(coord_df, exp_df) 
  }
}

st_get_assay = function(coord_df, exp_df, feature_col, feature_value_col, value_threshold ,top_n){
  	# Get feature and values
  	exp_df = exp_df[rownames(coord_df),]
    
    if(!missing(top_n)){
    	print(str_glue("selecting top {top_n} value from the {feature}"))
    	## Still thinking how to implement ....
    }
    if(!missing(value_threshold)){
    	print(str_glue('Filtering value by {value_threshold}'))
    	exp_df[exp_df < value_threshold] = 0 # set everything less than that threshold amount 0
    	exp_df[['alpha']] = apply(exp_df,1, sum, na.rm=T)
    }
    exp_df[is.na(exp_df)] = 0
    coord_df = cbind(coord_df[,(1:5)], exp_df)
}


filter_spots = function(coord_df, meta_df, filter_col, filter_group){
	# message(str_glue('nrow of coord_df = {nrow(coord_df)}'))
	selected_spot = meta_df %>% as.data.frame %>% rownames_to_column('barcode') %>%
		filter(!!as.symbol(all_of(filter_col)) %in% filter_group) %>% 
		pull('barcode')

	if(length(selected_spot) == 0){stop('No barcode left. Please check filtering parameter')}
	selected_spot = intersect(selected_spot, rownames(coord_df))
	coord_df = coord_df[selected_spot,]
	if(nrow(coord_df) == 0) message('No spot left. Please check filtering parameter')
	coord_df
}
###################################
######### Plot processing #########
st_scale_and_center = function(df, nodge_x=0, nodge_y=0){
  df %>% mutate(row = row * 127/78) %>%  # scale
    # mutate(row = row - 64+nodge_x, col = col-64+nodge_y) # center - v1(?)
    mutate(row = row - mean(row) , col = col - mean(col) ) # center - v2
}

st_rotation = function(coord_df, counter_degree = 90, mirror='none'){
  new_coord = coord_df[,c('row','col')] 
  # Rotation matrix (counter clock-wise)
  counter90  = matrix(c(0 , 1,-1, 0), ncol  = 2)
  counter180 = matrix(c(-1, 0, 0,-1), ncol  = 2)
  counter270 = matrix(c(0 ,-1, 1, 0),  ncol  = 2)
  
  new_coord = if(counter_degree==90) t(counter90  %*% t(new_coord)) 
  else if (counter_degree==180) t(counter180 %*% t(new_coord))
  else if (counter_degree==270) t(counter270 %*% t(new_coord))
  else new_coord # no rotation
  new_coord = new_coord %>% as.data.frame() %>% setNames(c('row','col'))
  coord_df$row = new_coord$row
  coord_df$col = new_coord$col

  if(mirror=='x') coord_df$row = -new_coord$row
  if(mirror=='y') coord_df$col = -new_coord$col
  if(mirror %in% c('xy','yx','both')) {
  	coord_df$row = -new_coord$row
  	coord_df$col = -new_coord$col
  }
  
  coord_df
}

## TESTING


st_ggplot = function(coord_df, feature, mode,
	xlimit = c(-64,64), ylimit = c(-64,64), dist_range, # distribution color range
	pt_size=1, xlab='', ylab='', show_axis =F, top_label='', meta_palette, alpha_col=NULL, gray_out=T,
	plot_spacing=0,
	boarder_size=0.25
	){

	# Highlight certain spots
	if(!is.null(alpha_col) & !gray_out){
		p = ggplot(coord_df, aes_string(x = 'row', y = 'col', fill = paste0("`",feature,"`"), alpha = alpha_col)) + 
			scale_alpha_continuous(limits = c(0,1))
	}else if(!is.null(alpha_col) & gray_out){ # Uses gray out mode
		# Split into 2 
		coord_gray_df = filter(coord_df, !!as.symbol(alpha_col) != 1) # unselected spots
		coord_df      = filter(coord_df, !!as.symbol(alpha_col) == 1) # highlight spots
		p = ggplot(coord_df, aes_string(x = 'row', y = 'col', fill = paste0("`",feature,"`")))
	}else{
		p = ggplot(coord_df, aes_string(x = 'row', y = 'col', fill = paste0("`",feature,"`")))
	}

  # Main plot 
  p = p + 
	    geom_point(size=pt_size, shape = 21, color = 'transparent') +
	    theme_void()+
	    labs(x = xlab, y = ylab) + 
	    coord_fixed() + 
	    NoLegend()
	
	# X/Y axis: Adjust spacing, limits and outline
	if(!show_axis) p = p + theme(axis.text=element_blank(), axis.line=element_blank(), axis.ticks=element_blank())

	p = p + 
			# Spacing in each plot
	    scale_x_continuous(expand=c(0,plot_spacing), limits = xlimit) + 
	    scale_y_continuous(expand=c(0,plot_spacing), limits = ylimit) +
	    # Add outline
	    theme(panel.border = element_rect(color = 'gray60', fill=NA, size=2),
	  				plot.margin  = unit(c(0,0,0,0), 'lines'))
	


	# Alternative rainbow color not using Seurat
	rainbow_palette <- colorRampPalette(rev(brewer.pal(11, "Spectral"))) # The same as Seurat:::SpatialColors

	# Plot mode switch 
	if(mode=='meta' & !missing(meta_palette)){ # Meta data plotting mode
		meta_palette = meta_palette[unique(coord_df[[feature]])]
		# message(meta_palette)
		message(str_glue('Meta mode using meta_palette of {length(meta_palette)}'))
		p = p + scale_fill_manual(values = meta_palette)

	# Feature plotting mode
	}else if(mode %in% c('proportion_adjust_feature','feature','meta_numeric', 'assay')){ 
		# No data 
		if(max(coord_df[[feature]])==0){
			#p = p + scale_fill_gradientn(colors = Seurat:::SpatialColors(n = 1)) 
			p = p + scale_fill_gradientn(colors = rainbow_palette(1)) #Seurat:::SpatialColors(n = 1)) 
		}else if(!is.null(dist_range)){
			#p = p + scale_fill_gradientn(colors = Seurat:::SpatialColors(n = 100), limits = dist_range) # fixed range
			p = p + scale_fill_gradientn(colors = rainbow_palette(100), limits = dist_range)# Seurat:::SpatialColors(n = 100), limits = dist_range) # fixed range
		}else{
			p = p + scale_fill_gradientn(colors = rainbow_palette(100)) #Seurat:::SpatialColors(n = 100)) 
		}
	}else{
		message(str_glue('Unexpected mode {mode}'))
	}

  # Add gray points if using gray out mode
	if(!is.null(alpha_col) & gray_out){
		p = p + geom_point(data = coord_gray_df, size = pt_size, shape = 21, fill = 'gray', color = 'transparent')
	}

	# Label and axis
	if(top_label!=''){
		p = p + 
			#scale_x_discrete(position = "top") +
			labs(title = top_label) + 
			theme(plot.title = element_text(size = 25, face = 'bold', angle =0, hjust =0.5))
	}
	p
}

# Highlight cells
st_highlight = function(df, obj, highlight_col, highlight_idents, highlight_alpha = 1, other_alpha = 0.1){
	highlight_df = obj@meta.data %>% mutate(highlights = ifelse(!!as.symbol(highlight_col) %in% highlight_idents, highlight_alpha, other_alpha))
	shared_spots = intersect(rownames(df), rownames(highlight_df))
	df[shared_spots,'highlights'] = highlight_df[shared_spots, 'highlights']
	df 
}

SingleSTDistIdentPlot = function(obj, slice, mode, meta_palette, rotation=0, feature, 
			pt_size=1, xlab='',ylab='', show_axis=F, 
			xlimit = c(-64,64), ylimit = c(-64,64),
			dist_range, # Distribution mode. Expression range from other plots
			assay_name, highlight_cell_type, proportion_adjust_feature = F,
			assay_use, # For plotting cell type proportion 'assay mode'
			nodge_x=0, nodge_y=0,mirror='none', top_label='', filter_meta, filter_col, filter_group, filter_piece_col, filter_piece, highlight_col, highlight_idents, gray_out=T){
	if(missing(mode)) stop("Need to specify either feature or meta mode for plotting")

	df = obj@images[[slice]]@coordinates

	# Highlight cells
	# V2 7/27/2021 - by mutiplying expression and cell proproption 
	if(proportion_adjust_feature){
		if(missing(highlight_cell_type)){
			stop("Using cell type proportion adjustment mode without cell type to highlight. Please use highlight_cell_type parameter to specify ")
		}
		message('Using proportion adjust feature mode')
		mode = 'proportion_adjust_feature'
	}
	# If highlight, make alpha_col
	# Need highlight col, highlight id
	if(!missing(highlight_col) & !missing(highlight_idents)){
		message(str_glue('Highlighting {toString(highlight_idents)} from {highlight_col}'))
		df = st_highlight(df, obj, highlight_col, highlight_idents, 1, 0.1)
		alpha_col = 'highlights'
	}else{# no highlights
		alpha_col = NULL
	}

	# Filter by gender or other variables 
	if(!missing(filter_meta)) df = df %>% filter_spots(meta_df =filter_meta, filter_col=filter_col, filter_group=filter_group) 
	if(missing(filter_piece_col)){filter_piece_col = NA; filter_piece=NA}
	if(!is.na(filter_piece_col)){
		df = df %>% filter_spots(meta_df =filter_meta, filter_col=filter_piece_col, filter_group=filter_piece)
		if(nrow(df)==0){
			avail_pieces = filter_meta %>% filter(orig.ident ==all_of(slice)) %>% pull('Piece_ID') %>% unique
			message(str_glue("Resulting df after filtering filter_piece is 0. \n {toString(filter_piece)} Piece_ID might be not in slice: {slice}"))
			message(str_glue("Available ones are {toString(avail_pieces)}"))
			stop('Please use m_filter_piece or f_filter_piece parameter to specify the right Piece_ID to filter E165')
			
		}
	}
	
	# distribution rnage
	if(missing(dist_range)) dist_range = NULL

	df = df %>% 
		  st_scale_and_center(nodge_x=nodge_x, nodge_y=nodge_y) %>%
		  st_rotation(rotation, mirror=mirror) 

	if(mode == 'feature'){
		df = df %>% st_get_expression(obj, feature) %>%
					st_ggplot(feature=feature, pt_size=pt_size, xlab=xlab, ylab=ylab, show_axis=show_axis, top_label = top_label, mode=mode, alpha_col = alpha_col, gray_out =gray_out, dist_range= dist_range) 

	}else if (mode =='proportion_adjust_feature'){ # 7/27/2021 - Adjust expression by selected cell type
		df = df %>% st_get_expression(obj, feature) %>%
						st_proportion_adjust_exp(assay_name, obj=st_merge, highlight_cell_type) %>%
						st_ggplot(feature=feature, pt_size=pt_size, xlab=xlab, ylab=ylab, show_axis=show_axis, top_label = top_label, mode=mode, alpha_col = alpha_col, gray_out =gray_out, dist_range= dist_range) 

	}else if (mode =='meta'){
		df = df %>% st_add_meta(obj, feature) 
		df = df %>% st_ggplot(feature=feature, pt_size=pt_size, xlab=xlab, ylab=ylab, show_axis=show_axis, top_label = top_label, mode=mode, meta_palette=meta_palette,
													xlimit = xlimit, ylimit = ylimit, alpha_col = alpha_col) 

	}else if(mode == 'meta_numeric'){
		df = df %>% st_add_meta(obj, feature) 
		df = df %>%
					st_ggplot(feature=feature, pt_size=pt_size, xlab=xlab, ylab=ylab, show_axis=show_axis, top_label = top_label, mode=mode, alpha_col = alpha_col, gray_out=gray_out, dist_range= dist_range) 
	}else if(mode =='assay'){
		message('assay mode hereeeee')
		df = df %>% st_get_assay_value(obj, feature, assay_use=assay_use) %>%
					st_ggplot(feature=feature, pt_size=pt_size, xlab=xlab, ylab=ylab, show_axis=show_axis, top_label = top_label, mode=mode, alpha_col = alpha_col, gray_out =gray_out, dist_range= dist_range) 
	}
	else{
		message('Only support mode = feature or meta')
	}
	df
}

# 7/27/2021: Adjust expression to highlihgt cell by multiple cell type proportion with expression level
st_proportion_adjust_exp = function(coord_exp_df, obj, assay_name = 'RCTD_doubletRM_v3', highlight_cell_type){
	if(!assay_name %in% Assays(obj)){
		stop(str_glue('Assay {assay_name} not found in the object. Available ones are {toString(Assays(st_merge))}'))
	}
	if(missing(highlight_cell_type)){
		stop("Missing cell type to highlight. Please make sure to specify")
	}
	avail_cell_types = rownames(GetAssay(obj, assay = assay_name) %>% GetAssayData())
	if(!highlight_cell_type %in% avail_cell_types){
		stop(str_glue('{highlight_cell_type} not found in assay {assay_name}. Available ones are {toString(avail_cell_types)}'))
	}
	message(str_glue('Adjusting expresion to {highlight_cell_type} proportion'))
	# Get Cell type porportion
	celltype_proportions = GetAssay(obj, assay = assay_name) %>% GetAssayData()
	celltype_proportions = celltype_proportions %>% .[highlight_cell_type, rownames(coord_exp_df)] 
	# Check if all NA
	if(all(is.na(celltype_proportions))){
		message(str_glue('All {highlight_cell_type} in {assay_name} is NA, aka no cell type detected in this slice. Use 0 for all values.'))
		celltype_proportions = 0
	}
	
	# Get gene name
	gene = setdiff(colnames(coord_exp_df), c("tissue", "row", "col", "imagerow", "imagecol"))
	
	# Adjust expression by multiply expression with ratio 
	#coord_exp_df[,str_c(gene,'_orig')] = coord_exp_df[,gene] # Keep original for comparison
	coord_exp_df[,gene] = coord_exp_df[,gene] * celltype_proportions
	coord_exp_df
}


MultipleSTDistIdentPlot = function(obj, images, feature, mode, left_label = '', mirrors, rotates, top_labels, nodges_x, nodges_y, filter_piece_col , filter_piece, xlimits, ylimits, ...){
	if(missing(mirrors)) mirrors = rep('no', length(images)) %>% setNames(images)
	if(missing(rotates)) rotates = rep(0, length(images)) %>% setNames(images)
	if(missing(top_labels)) top_labels = rep("", length(images)) %>% setNames(images)
	if(missing(nodges_x)) nodges_x = rep(0, length(images)) %>% setNames(images)
	if(missing(nodges_y)) nodges_y = rep(0, length(images)) %>% setNames(images)
	if(missing(filter_piece_col))   filter_piece_col    = NA # Use only in the E165 for the second filtering/spot selection
	if(missing(filter_piece)) filter_piece  = NA # Use only in the E165 for the second filtering/spot selection
	if(missing(xlimits)) xlimits = rep(list(c(-64,64)), length(images)) %>% setNames(images) # Named list with x limits
	if(missing(ylimits)) ylimits = rep(list(c(-64,64)), length(images)) %>% setNames(images) # Named list with y limits
	


	p1     = SingleSTDistIdentPlot(obj=obj, slice=images[[1]], feature = feature, mode = mode, 
								  mirror    = mirrors[[1]],    rotation = rotates[[1]], 
								  top_label = top_labels[[1]], nodge_x=nodges_x[[1]], nodge_y=nodges_y[[1]], 
								  filter_piece_col =filter_piece_col, filter_piece =filter_piece, 
								  xlimit = xlimits[[images[[1]]]], ylimit = ylimits[[images[[1]]]], 
								  ...) + 
				labs(y=left_label) + 
				theme(axis.title.y = element_text(size = 25, face = 'bold', angle =0, vjust =0.5))

	p_rest = map(images[-1], function(slice){
		SingleSTDistIdentPlot(obj=obj, slice=slice, feature = feature, mode = mode, 
							  mirror = mirrors[[slice]], rotation = rotates[[slice]], 
							  top_label = top_labels[[slice]], nodge_x=nodges_x[[slice]], nodge_y=nodges_y[[slice]],
							  filter_piece_col = NA, filter_piece = NA,
							  xlimit = xlimits[[slice]], ylimit = ylimits[[slice]], 
							  ...)
	})
	c(list(p1), p_rest)
}


cell_type_scatterpie = function(obj, slice, cell_type_df, value_threshold=0.1, rotation=0, pie_size=1, xlimit=c(-64,64), ylimit = c(-64,64), palette, use_alpha=T){
	# Make df
	pie_df = obj@images[[slice]]@coordinates %>% 
	  st_scale_and_center %>%
	  st_get_assay(st_all_data, value_threshold=value_threshold)  %>% 
	  st_rotation(rotation) %>% 
	  mutate(radius = pie_size) 
	
	# Select cols
	cols_to_plot = setdiff(names(pie_df), c("tissue", "row", "col", "imagerow", "imagecol", "alpha", "radius"))
	
	# Make pie plot
	p = ggplot() 
	if(use_alpha) p = p + geom_scatterpie(data = pie_df, aes(x = row, y = col, r=radius, alpha = alpha), cols = cols_to_plot,
                     color = NA)
	else p = p + geom_scatterpie(data = pie_df, aes(x = row, y = col, r=radius), cols = cols_to_plot,
                     color = NA)
    
    p = p + xlim(xlimit) +
	    ylim(ylimit) +
	    coord_fixed() + 
	    theme_bw() + 
	    theme(panel.grid = element_blank())

    if(!missing(palette)) p = p + scale_fill_manual(values = palette)

    p
}
### EXAMPLE ####
# Cell type Pie plot with filtering
# pdf('figure/manuscript/fig1/Cell_Type_filter_pie_01.pdf')
# cell_type_scatterpie(st_merge, 'NMK0205E165MU2', st_all_data, palette = col_celltype)
# dev.off()


