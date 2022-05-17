## ST NMK plotting function ###
## 06/09/2021 ####
## Simon Mo ######
source('/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/NMK/NMK_integration_v5/script/NMK_colorset_v3_07152021.R')
source('/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/NMK/NMK_integration_v5/script/Functions/function_ST_plot_general_v1.R')


make_nmk_mf_st_plot = function(...){
	message('make_nmk_mf_st_plot function is being deprecated. Use STDistPlot_NMK_MF instead')
	STDistPlot_NMK_MF(...)
}

###### ST Plotting functions -11022021 #####
# internal
STDistAKIPlots = function(obj, features, ...){
    map(features,  function(gene) STDistPlot_AKI(obj, feature = gene, show_axis=F,...)) %>% setNames(features)
}
STDistNMKPlots = function(obj, features, ...){
    map(features,  function(gene) STDistPlot_NMK_MF(obj, 
                                      feature = gene, 
                                      m_filter_piece = c('NMK0205E165M2U1L','NMK0205E165M2U1R'),
                                      ...)) %>% setNames(features)
    }

# API
STDistPlots = function(obj, features, mode = c('AKI','NMK'), pt_size=1.8, ...){
    mode   = match.arg(mode)
    pltfxn = switch(mode, 
                    'AKI'=STDistAKIPlots, 
                    'NMK'=STDistNMKPlots)
    pltfxn(obj, features, pt_size=pt_size, ...)

    }
STDistPlotsBoth = function(st_nmk, st_aki, features, pt_size=1.8, combine=T, ...){
    pNMKs = STDistPlots(st_nmk, features, mode = 'NMK', pt_size=pt_size, ...) %>% setNames(features)
    pAKIs = STDistPlots(st_aki, features, mode = 'AKI', pt_size=pt_size, ...) %>% setNames(features)
        if(combine){
        map(features, function(gene) 
            wrap_plots(N=pNMKs[[gene]], A=pAKIs[[gene]], design = "NNNNNNNNAAA") + 
                plot_annotation(title = gene, theme = theme(plot.title = element_text(size=20, face='bold')))
            ) %>% setNames(features)
        }else{
            list('NMK'=pNMKs, 'AKI'=pAKIs)
        }
    }
getSTPlotSize = function(type=c('NMK','AKI','Both')){
    type = match.arg(type) 
    switch(type, 
           'NMK' =c(35,10), 
           'AKI' =c(9,9), 
           'Both'=c(40,10)) # plot size
    }
PlotSTJupyter = function(p, type=c('NMK','AKI','Both')){
    type = match.arg(type) 
    pdim = getSTPlotSize(type)
    
    options(repr.plot.width=pdim[[1]], repr.plot.height=pdim[[2]])
    suppressWarnings(plot(p))
    }
PlotSTlistJupyter = function(plist, type){
    walk(plist, function(p) PlotSTJupyter(p, type))
}
SaveSTPlots = function(plist, type=c('NMK','AKI','Both'), path){
    today_date = Sys.Date() %>% format('%m%d%Y')
    type = match.arg(type) 
    pdim = getSTPlotSize(type)
    
    iwalk(plist, function(p, gene){
        ggsave(plot=p, str_glue("{path}/STDist_{gene}_{today_date}.pdf"),
               width = pdim[[1]],height = pdim[[2]]) %>% suppressWarnings
    })
    print(str_glue('Finished output {length(plist)} plots to {path}'))
}
#######
         

## STExpression value range 
getSTExpRange = function(obj, imgs, feature, mode, assay_use){
	cells      = map(imgs, function(img){obj@images[[img]]@coordinates %>% rownames}) %>% unlist 
	if(mode %in% c('feature')){
		exp_range = range( GetAssayData(obj) %>% .[feature, cells, drop=F] ) 
	}else if(mode =='meta_numeric'){
		exp_range = range( obj@meta.data[cells, feature, drop=F])
	}else if(mode =='assay'){
		exp_range = GetAssay(obj, assay = assay_use) %>% GetAssayData() %>% t %>% .[cells, feature, drop=F] %>% range
	}else{
		stop("Mode not supported. Please double check")
	}
	exp_range
}

## AKI version 10062021
STDistPlot_AKI = function(obj, feature, pt_size=1.6, images_m, images_f, xlimits_m, ylimits_m, xlimits_f, ylimits_f, 
		meta_palette,
		assay_use='RCTD_doubletRM_v3',
		plot_title, plot_title_size=30, plot_subtitle_size=20, plot_caption_size=15,
		plot_subtitle,
		plot_caption,
		exp_range_mode = 'all',
		...
	){
	if(!missing(images_m)){
		print('Warning. Using custom image from AKI object. All other parameters might not work well.')
	}else{
		images_m  = c('NMK0129M3MU1','NMK0129AKID3MU1','NMK0203AKID8MU1')
	}
	if(!missing(images_f)){
		print('Warning. Using custom image from AKI object. All other parameters might not work well.')
	}else{
		images_f  = c('NMK0129M3FU1','NMK0129AKID3FU1','NMK1225AKID8FU1')
	}

	# Choose mode
	if(feature %in% rownames(obj)){
		message(str_glue('{feature} found in expression table. Use feature mode.')); mode = 'feature'
	}else if(feature %in% colnames(obj@meta.data)){
		if(is.numeric(obj@meta.data[[feature]])){
			message(str_glue('{feature} found in @meta.data but numeric. Use meta_numeric mode.')); mode = 'meta_numeric'
		}else{
			message(str_glue('{feature} found in @meta.data. Use meta mode.')); mode = 'meta'
		}
	}else if(feature %in% rownames(GetAssay(obj, assay = assay_use))){
		message(str_glue('{feature} found in assay {assay_use}. Use assay mode')); mode = 'assay'
	}else{
		stop(str_glue('{feature} not found in expression, meta or assay {assay_use}. please double check.'))
	}

	# Check if have palette for meta mode
	if(mode =='meta' & missing(meta_palette)){
		message("Plotting meta but missing meta_palette. Use defualt color")
		meta_palette = pals::cols25(n = length(unique(obj@meta.data[[feature]]))) %>% setNames(unique(obj@meta.data[[feature]]))
	}

	# Set parameter
	if(missing(xlimits_m)) xlimits_m = rep(list(c(-64,64)), length(images_m)) %>% setNames(images_m) # Named list with x limits
	if(missing(xlimits_f)) xlimits_f = rep(list(c(-64,64)), length(images_f)) %>% setNames(images_f) # Named list with x limits
	if(missing(ylimits_m)) ylimits_m = rep(list(c(-64,64)), length(images_m)) %>% setNames(images_m) # Named list with x limits
	if(missing(ylimits_f)) ylimits_f = rep(list(c(-64,64)), length(images_f)) %>% setNames(images_f) # Named list with x limits

	## Parameters v1 
	## 06/09/2021
	
	label_m   = c('Control','AKI-Day3','AKI-Day8') %>% setNames(images_m)
	label_f   = c('Control','AKI-Day3','AKI-Day8') %>% setNames(images_f)
	# label_f   = rep('',8) %>% setNames(images_f)
	rotate_m  = c(0,0,0) %>% setNames(images_m)
	rotate_f  = c(0,0,0) %>% setNames(images_f)
	nodge_x_m = c(0,0,0) %>% setNames(images_m)
	nodge_y_m = c(0,0,0) %>% setNames(images_m)
	nodge_x_f = c(0,0,0) %>% setNames(images_f)
	nodge_y_f = c(0,0,0) %>% setNames(images_f)
	mirror_m  = c('no','no','no') %>% setNames(images_m)
	mirror_f  = c('no','no','no') %>% setNames(images_f)

	## Expresssion Range mode
	if(exp_range_mode == 'all'){
		message("Exp color range set to All sample")
		dist_range_f = getSTExpRange(obj,  c(images_m, images_f), feature=feature, mode = mode)
		dist_range_m = getSTExpRange(obj, c(images_m, images_f), feature=feature, mode = mode)
	}else if(exp_range_mode == 'by_sex'){
		message("Exp color range set by Sex")
		dist_range_f = getSTExpRange(obj, c(images_f), feature=feature, mode = mode)
		dist_range_m = getSTExpRange(obj, c(images_m), feature=feature, mode = mode)
	}else{
		message("No Exp color range set")
		dist_range_f = NULL
		dist_range_m = NULL
	}

	


	#### Plot gene features #####
	if(mode == 'feature'){
		p_mlist = MultipleSTDistIdentPlot(obj, images_m, feature= feature, left_label='M', mode = mode,
										pt_size=pt_size, mirror =mirror_m, top_labels=label_m,
										filter_meta = obj@meta.data, filter_col='Sex', 
										filter_group='M', 
										nodges_x= nodge_x_m, nodges_y= nodge_y_m,
										xlimits = xlimits_m, ylimits = ylimits_m,
										dist_range = dist_range_m, # Pass to SingleSTDistIdentPlot. Expression range
										...
										)
		p_flist = MultipleSTDistIdentPlot(obj, images_f, feature= feature, left_label='F', mode = mode,
										pt_size=pt_size, mirror =mirror_f, #top_labels=label_f,
										filter_meta = obj@meta.data, filter_col='Sex', 
										filter_group='F',
										nodges_x= nodge_x_f, nodges_y= nodge_y_f,
										xlimits = xlimits_f, ylimits = ylimits_f,
										dist_range = dist_range_f, # Pass to SingleSTDistIdentPlot. Expression range
										...
										)
	}else if(mode =='meta'){
		p_mlist = MultipleSTDistIdentPlot(obj, images_m, feature= feature, left_label='M', mode = mode, meta_palette = meta_palette,
										pt_size=pt_size, mirror =mirror_m, top_labels=label_m,
										filter_meta = obj@meta.data, filter_col='Sex', 
										filter_group='M', 
										nodges_x= nodge_x_m, nodges_y= nodge_y_m,
										xlimits = xlimits_m, ylimits = ylimits_m,
										...
										)
		p_flist = MultipleSTDistIdentPlot(obj, images_f, feature= feature, left_label='F', mode = mode, meta_palette = meta_palette,
										pt_size=pt_size, mirror =mirror_f, #top_labels=label_f,
										filter_meta = obj@meta.data, filter_col='Sex', 
										filter_group='F',
										nodges_x= nodge_x_f, nodges_y= nodge_y_f,
										xlimits = xlimits_f, ylimits = ylimits_f,
										...
										)
	}else if(mode =='meta_numeric'){
		p_mlist = MultipleSTDistIdentPlot(obj, images_m, feature= feature, left_label='M', mode = mode,
										pt_size=pt_size, mirror =mirror_m, top_labels=label_m,
										filter_meta = obj@meta.data, filter_col='Sex', 
										filter_group='M', 
										nodges_x= nodge_x_m, nodges_y= nodge_y_m,
										xlimits = xlimits_m, ylimits = ylimits_m,
										dist_range = dist_range_m, # Pass to SingleSTDistIdentPlot. Expression range
										...
										)
		p_flist = MultipleSTDistIdentPlot(obj, images_f, feature= feature, left_label='F', mode = mode,
										pt_size=pt_size, mirror =mirror_f, #top_labels=label_f,
										filter_meta = obj@meta.data, filter_col='Sex', 
										filter_group='F',
										nodges_x= nodge_x_f, nodges_y= nodge_y_f,
										xlimits = xlimits_f, ylimits = ylimits_f,
										dist_range = dist_range_f, # Pass to SingleSTDistIdentPlot. Expression range
										...
										)
	}else if(mode == 'assay'){ # 7/27/2021 For cell assay
		p_mlist = MultipleSTDistIdentPlot(obj, images_m, feature= feature, left_label='M', mode = mode,
										pt_size=pt_size, mirror =mirror_m, top_labels=label_m,
										filter_meta = obj@meta.data, filter_col='Sex', 
										filter_group='M', 
										nodges_x= nodge_x_m, nodges_y= nodge_y_m,
										xlimits = xlimits_m, ylimits = ylimits_m,
										assay_use = assay_use,
										...
										)
		p_flist = MultipleSTDistIdentPlot(obj, images_f, feature= feature, left_label='F', mode = mode,
										pt_size=pt_size, mirror =mirror_f, #top_labels=label_f,
										filter_meta = obj@meta.data, filter_col='Sex', 
										filter_group='F',
										nodges_x= nodge_x_f, nodges_y= nodge_y_f,
										xlimits = xlimits_f, ylimits = ylimits_f,
										assay_use = assay_use,
										...
										)
	}

	# Combine plots
	male_widths = map(xlimits_m, ~abs(.x[[1]]-.x[[2]])) %>% unlist
	male_hights = map(ylimits_m, ~abs(.x[[1]]-.x[[2]])) %>% unlist
	switch(missing(plot_subtitle), NULL, plot_subtitle)
	message('Current version works well with width = ?, height = 10 ')
	mstrip = wrap_plots(p_mlist, nrow=1, widths = male_widths, heights =male_hights)
	fstrip = wrap_plots(p_flist, nrow=1, widths = male_widths, heights =male_hights)

	# Plot annotations
	if(missing(plot_title))      plot_title = feature
	if(missing(plot_subtitle))   plot_subtitle = ""
	if(missing(plot_caption))    plot_caption = ""

	((mstrip/fstrip) &
	theme(plot.margin = unit(c(0,0,0,0),'lines')) )+ ## Still dont know how to remove space between 2 stirps ....
	plot_annotation(
		title = plot_title,
		subtitle = plot_subtitle,
		caption = plot_caption,
		theme = theme(plot.title = element_text(face = 'bold', size = plot_title_size),
					  plot.subtitle = element_text(face = 'bold', size = plot_subtitle_size),
					  plot.caption = element_text(face = 'bold', size = plot_caption_size))
	)
}



STDistPlot_NMK_MF = function(obj, feature, pt_size=1.6, images_m, images_f, xlimits_m, ylimits_m, xlimits_f, ylimits_f, 
		m_filter_piece = c('NMK0205E165M2U1L','NMK0205E165M2U1R'), # Use on the first image for m
		f_filter_piece = c('NMK0318E165F2U1L','NMK0318E165F2U1R'),  # Use on the first image for m,
		m_filter_group1, f_filter_group1,
		meta_palette,
		assay_use='RCTD_doubletRM_v3',
		plot_title, plot_title_size=30, plot_subtitle_size=20, plot_caption_size=15,
		plot_subtitle,
		plot_caption,
		...
	){
	if(!missing(images_m)){
		print('Warning. Using custom image from NMK object. All other parameters might not work well.')
	}else{
		images_m  = c('NMK0205E165MU1','NMK0117P0FMU1L','NMK0115W1FMU1','NMK0111W2MU1','NMK1201W3MU1','NMK0129W12MU1','NMK1207W52MU1','NMK0226W113MU1')
	}
	if(!missing(images_f)){
		print('Warning. Using custom image from NMK object. All other parameters might not work well.')
	}else{
		images_f  = c('NMK0318E165FU1','NMK0117P0FMU1L','NMK0115W1FMU1','NMK0111W2FU1','NMK1201W3FU2','NMK0129W12FU1','NMK1207W52FU1','NMK0430W113FU1')
	}

	# Choose mode
	if(feature %in% rownames(obj)){
		message(str_glue('{feature} found in expression table. Use feature mode.')); mode = 'feature'
	}else if(feature %in% colnames(obj@meta.data)){
		if(is.numeric(obj@meta.data[[feature]])){
			message(str_glue('{feature} found in @meta.data but numeric. Use meta_numeric mode.')); mode = 'meta_numeric'
		}else{
			message(str_glue('{feature} found in @meta.data. Use meta mode.')); mode = 'meta'
		}
	}else if(feature %in% rownames(GetAssay(obj, assay = assay_use))){
		message(str_glue('{feature} found in assay {assay_use}. Use assay mode')); mode = 'assay'
	}else{
		stop(str_glue('{feature} not found in expression, meta or assay {assay_use}. please double check.'))
	}

	# Check if have palette for meta mode
	if(mode =='meta' & missing(meta_palette)){
		message("Plotting meta but missing meta_palette. Use defualt color")
        meta_palette = paletteer::paletteer_d("ggsci::default_igv", n = length(unique(obj@meta.data[[feature]]))) %>% setNames(unique(obj@meta.data[[feature]]))
		# meta_palette = pals::cols25(n = length(unique(obj@meta.data[[feature]]))) %>% setNames(unique(obj@meta.data[[feature]]))
	}

	# xlimts Parameters v06220221
	#if(missing(xlimits_m)) xlimits_m = c(list(c(-40,40)), list(c(-22,21)) ,list(c(-35,30)), list(c(-55,64)), rep(list(c(-64,64)), 8-4)) %>% setNames(images_m)
	#if(missing(xlimits_f)) xlimits_f = c(list(c(-40,40)), list(c(-22,21)) ,list(c(-35,30)), list(c(-55,64)), rep(list(c(-64,64)), 8-4)) %>% setNames(images_f)


	# !!! DEPRECATING filter group become filter piece
	if(!missing(m_filter_group1)) {
		warning("m_filter_group1 being deprecated. Use m_filter_piece instead")
		m_filter_piece = m_filter_group1
	}
	if(!missing(f_filter_group1)) {
		warning("f_filter_group1 being deprecated. Use f_filter_piece instead")
		f_filter_piece = f_filter_group1
	}

	# Set parameter
	if(missing(xlimits_m)) xlimits_m = rep(list(c(-64,64)), length(images_m)) %>% setNames(images_m) # Named list with x limits
	if(missing(xlimits_f)) xlimits_f = rep(list(c(-64,64)), length(images_f)) %>% setNames(images_f) # Named list with x limits
	if(missing(ylimits_m)) ylimits_m = rep(list(c(-64,64)), length(images_m)) %>% setNames(images_m) # Named list with x limits
	if(missing(ylimits_f)) ylimits_f = rep(list(c(-64,64)), length(images_f)) %>% setNames(images_f) # Named list with x limits

	## Parameters v1 
	## 06/09/2021
	
	label_m   = c('E165','P0','W1','W2','W3','W12','W52','W113') %>% setNames(images_m)
	label_f   = c('E165','P0','W1','W2','W3','W12','W52','W113') %>% setNames(images_f)
	# label_f   = rep('',8) %>% setNames(images_f)
	rotate_m  = c(0,0,0,0,0,0,0,0) %>% setNames(images_m)
	rotate_f  = c(0,0,0,0,0,0,0,0) %>% setNames(images_f)
	nodge_x_m = c(0,0,0,0,0,0,0,0) %>% setNames(images_m)
	nodge_y_m = c(-20,10,0,10,0,0,0,0) %>% setNames(images_m)
	nodge_x_f = c(0,0,0,0,0,0,0,0) %>% setNames(images_f)
	nodge_y_f = c(-35,5,-10,0,0,0,0,0) %>% setNames(images_f)
	mirror_m  = c('no','xy','x','no','no','no','y','no') %>% setNames(images_m)
	mirror_f  = c('no','no','no','no','no','no','y','xy') %>% setNames(images_f)


	## Run these to check and select E165 pieces in 'filter_piece'
	#obj@meta.data %>% filter(orig.ident ==images_m[[1]]) %>% count(Piece_ID)
	#obj@meta.data %>% filter(orig.ident ==images_f[[1]]) %>% count(Piece_ID)
	##

	## Get distribution value range 
	cells      = map(c(images_m, images_f), function(img){obj@images[[img]]@coordinates %>% rownames}) %>% unlist 
	if(mode %in% c('feature')){
		dist_range = range( GetAssayData(obj) %>% .[feature, cells, drop=F] ) 
	}else if(mode =='meta_numeric'){
		dist_range = range( obj@meta.data[cells, feature, drop=F])
	}


	#### Plot gene features #####
	if(mode == 'feature'){
		p_mlist = MultipleSTDistIdentPlot(obj, images_m, feature= feature, left_label='M', mode = mode,
										pt_size=pt_size, mirror =mirror_m, top_labels=label_m,
										filter_meta = obj@meta.data, filter_col='Sex', 
										filter_group='M', 
										nodges_x= nodge_x_m, nodges_y= nodge_y_m,
										filter_piece_col = 'Piece_ID', filter_piece = m_filter_piece,
										xlimits = xlimits_m, ylimits = ylimits_m,
										dist_range = dist_range, # Pass to SingleSTDistIdentPlot
										...
										)
		p_flist = MultipleSTDistIdentPlot(obj, images_f, feature= feature, left_label='F', mode = mode,
										pt_size=pt_size, mirror =mirror_f, #top_labels=label_f,
										filter_meta = obj@meta.data, filter_col='Sex', 
										filter_group='F',
										nodges_x= nodge_x_f, nodges_y= nodge_y_f,
										filter_piece_col = 'Piece_ID', filter_piece = f_filter_piece,
										xlimits = xlimits_f, ylimits = ylimits_f,
										dist_range = dist_range, # Pass to SingleSTDistIdentPlot
										...
										)
	}else if(mode =='meta'){
		p_mlist = MultipleSTDistIdentPlot(obj, images_m, feature= feature, left_label='M', mode = mode, meta_palette = meta_palette,
										pt_size=pt_size, mirror =mirror_m, top_labels=label_m,
										filter_meta = obj@meta.data, filter_col='Sex', 
										filter_group='M', 
										nodges_x= nodge_x_m, nodges_y= nodge_y_m,
										filter_piece_col = 'Piece_ID', filter_piece = m_filter_piece,
										xlimits = xlimits_m, ylimits = ylimits_m,
										...
										)
		p_flist = MultipleSTDistIdentPlot(obj, images_f, feature= feature, left_label='F', mode = mode, meta_palette = meta_palette,
										pt_size=pt_size, mirror =mirror_f, #top_labels=label_f,
										filter_meta = obj@meta.data, filter_col='Sex', 
										filter_group='F',
										nodges_x= nodge_x_f, nodges_y= nodge_y_f,
										filter_piece_col = 'Piece_ID', filter_piece = f_filter_piece,
										xlimits = xlimits_f, ylimits = ylimits_f,
										...
										)
	}else if(mode =='meta_numeric'){
		p_mlist = MultipleSTDistIdentPlot(obj, images_m, feature= feature, left_label='M', mode = mode,
										pt_size=pt_size, mirror =mirror_m, top_labels=label_m,
										filter_meta = obj@meta.data, filter_col='Sex', 
										filter_group='M', 
										nodges_x= nodge_x_m, nodges_y= nodge_y_m,
										filter_piece_col = 'Piece_ID', filter_piece = m_filter_piece,
										xlimits = xlimits_m, ylimits = ylimits_m,
										...
										)
		p_flist = MultipleSTDistIdentPlot(obj, images_f, feature= feature, left_label='F', mode = mode,
										pt_size=pt_size, mirror =mirror_f, #top_labels=label_f,
										filter_meta = obj@meta.data, filter_col='Sex', 
										filter_group='F',
										nodges_x= nodge_x_f, nodges_y= nodge_y_f,
										filter_piece_col = 'Piece_ID', filter_piece = f_filter_piece,
										xlimits = xlimits_f, ylimits = ylimits_f,
										...
										)
	}else if(mode == 'assay'){ # 7/27/2021 For cell assay
		p_mlist = MultipleSTDistIdentPlot(obj, images_m, feature= feature, left_label='M', mode = mode,
										pt_size=pt_size, mirror =mirror_m, top_labels=label_m,
										filter_meta = obj@meta.data, filter_col='Sex', 
										filter_group='M', 
										nodges_x= nodge_x_m, nodges_y= nodge_y_m,
										filter_piece_col = 'Piece_ID', filter_piece = m_filter_piece,
										xlimits = xlimits_m, ylimits = ylimits_m,
										assay_use = assay_use,
										...
										)
		p_flist = MultipleSTDistIdentPlot(obj, images_f, feature= feature, left_label='F', mode = mode,
										pt_size=pt_size, mirror =mirror_f, #top_labels=label_f,
										filter_meta = obj@meta.data, filter_col='Sex', 
										filter_group='F',
										nodges_x= nodge_x_f, nodges_y= nodge_y_f,
										filter_piece_col = 'Piece_ID', filter_piece = f_filter_piece,
										xlimits = xlimits_f, ylimits = ylimits_f,
										assay_use = assay_use,
										...
										)
	}

	# Combine plots
	male_widths = map(xlimits_m, ~abs(.x[[1]]-.x[[2]])) %>% unlist
	male_hights = map(ylimits_m, ~abs(.x[[1]]-.x[[2]])) %>% unlist
	switch(missing(plot_subtitle), NULL, plot_subtitle)
	message('Current version works well with width = 35, height = 10 ')
	mstrip = wrap_plots(p_mlist, nrow=1, widths = male_widths, heights =male_hights)
	fstrip = wrap_plots(p_flist, nrow=1, widths = male_widths, heights =male_hights)

	# Plot annotations
	if(missing(plot_title))      plot_title = feature
	if(missing(plot_subtitle))   plot_subtitle = ""
	if(missing(plot_caption))    plot_caption = ""

	((mstrip/fstrip) &
	theme(plot.margin = unit(c(0,0,0,0),'lines')) )+ ## Still dont know how to remove space between 2 stirps ....
	plot_annotation(
		title = plot_title,
		subtitle = plot_subtitle,
		caption = plot_caption,
		theme = theme(plot.title = element_text(face = 'bold', size = plot_title_size),
					  plot.subtitle = element_text(face = 'bold', size = plot_subtitle_size),
					  plot.caption = element_text(face = 'bold', size = plot_caption_size))
	)
}



make_nmk_multiple_st_plot = function(obj, features, ...){
	p_list = map(features, function(feature){
		STDistPlot_NMK_MF(obj = obj, feature = feature, ...)
	}) %>% setNames(features)
	message("Done. Output as a list of plots")
	p_list
}
