# Spatial plotting helper
library(Seurat)
library(tidyverse)
library(patchwork)
library(viridis)

# Helper function
# Use image name to replace orig.ident for plotting
ReplaceMetaWithImageName = function(obj, target_col = 'orig.ident'){
  # Get Image Name Table
  image_name_table = map(Images(obj), function(slide_name){
    data.frame(id = obj@images[[slide_name]]@coordinates %>% rownames(), image_name = slide_name)
  }) %>% bind_rows() %>% column_to_rownames(('id'))
  # Replace orig.ident in obj
  AddMetaData(obj, metadata = image_name_table[rownames(obj@meta.data), 'image_name', drop=T], col.name = target_col)
}



SpatialDimPlotOrganizedList = function(obj, plot_meta,  nrow = 3, rainbow_col =T, hide_image=F){

map(plot_meta, function(group.by){
    if(hide_image){
        p = SpatialPlot(obj, group.by = group.by, combine =F, stroke=NA, image.alpha = 0) 
    }else{
        p = SpatialPlot(obj, group.by = group.by, combine =F, stroke=NA) 
    }
	# Main
	p = p %>% wrap_plots(guides = 'collect', nrow=nrow) & 
    	theme(legend.box = 'horizontal', legend.position = "bottom") &
    	guides(fill = guide_legend(ncol=2)) 

    # Color
    if(rainbow_col){
		color_map = get_rainbow_col(length(unique(obj@meta.data[[group.by]])) ) %>% 
					              setNames(unique(obj@meta.data[[group.by]]) )
		p = p & scale_fill_manual(values = color_map)
	}

    # Annotation
    p + plot_annotation(title = group.by)
})

}



SpatialFeaturePlotPlus = function(obj, features, images, expression_ref_obj, sample_labels, pt.size.factors, feature_position = 'column', image.alpha = 0, stroke=NA, color_palette = 'SpatialColor', assemble=T, ...){
    # Images/sample to plot
    if(missing(images)){
        message(str_glue("Using all {length(Images(obj))} images available in the Spatial Seurat obj"))
        images = Images(obj)
    }
    if(!all(images %in% obj$orig.ident)){
        # Try fix it
        obj = ReplaceMetaWithImageName(obj)
        if(!all(images %in% obj$orig.ident)){ # check again
            stop("Not all images provided found in orig.ident. Make sure image names are the same in orig.ident")
        }
    }
    # # Subset object with only given sample and features
    
    # Idents(obj) = 'orig.ident'
    # if(!identical( sort(images), sort( unique(obj$orig.ident) ) )){ # Not using all sample - Subset
    #     message('Subsetting obj..')
    #     obj = subset(obj, orig.ident %in% images)[features, ]
    # }
    if(!all(features%in%rownames(obj))){
        stop(str_glue("Error: {toString(setdiff(features, rownames(obj)))} features not found in the assay"))
    }
    # features  = intersect( features, rownames(obj))

    # Handle color scale limit
    message('Getting lower and upper bound of each feature..')
    if(!missing(expression_ref_obj)){
        message('Getting expression range from a separated expression obj')
        exp_df = GetAssayData(expression_ref_obj)
    }else{
        message('Getting expression range using the same obj for plotting')
        exp_df = GetAssayData(obj)
    }
    exp_df    = exp_df[features,]
    if(is.null(dim(exp_df))){ # only 1 feature
        exp_lows  = min(exp_df)
        exp_highs = max(exp_df)
    }else{
        exp_lows  = exp_df %>% apply(1, min)
        exp_highs = exp_df %>% apply(1, max)
    }
    
    # Handle missing vectors
    # Handle sample label on the left - for sample name on the left of plot
    if(missing(sample_labels)){
        sample_labels = images
    }else if(length(sample_labels)!=length(images)){
        message('Warning: sample_labels lenght not equal to images length. Label with images names')
        sample_labels = images
    }

    # pt size
    if(missing(pt.size.factors)){
        defualt_ptsize  = 1.6 
        pt.size.factors = rep(defualt_ptsize, length(images))
    }else if(length(pt.size.factors)!=length(images)){
        message('Warning: pt.size.factors lenght not equal to images length. Using first element only')
        pt.size.factors = rep(pt.size.factors[[1]], length(images))
    }

    

    # Choose color palette
    col_palette = case_when(
        color_palette=='viridis' ~ viridis(n=100),
        color_palette=='magma' ~ magma(n=100),
        color_palette=='plasma' ~ plasma(n=100),
        color_palette=='inferno' ~ inferno(n=100),
                            TRUE ~ Seurat:::SpatialColors(n = 100)
        )

    # Choose plot type
    if(feature_position == 'column'){
        message("Feature as column, sample as row.")
        p = SpatialFeaturePlotPlusSampleAsRow(obj, features, images, sample_labels, pt.size.factors, exp_lows, exp_highs, col_palette, image.alpha = 0, stroke=NA, assemble_col =assemble, ...)
    }
    else{
        message("Feature as row, sample as column.")
        p = SpatialFeaturePlotPlusFeatureAsRow(obj, features, images, sample_labels, pt.size.factors, exp_lows, exp_highs, col_palette, image.alpha = 0, stroke=NA, assemble_row =assemble,  ...)
    }
    p

}

SpatialFeaturePlotPlusSampleAsRow = function(obj, features, images, sample_labels, pt.size.factors, exp_lows, exp_highs, col_palette, assemble_col=T, ...){
    # Seurat plotting function
    message("Creating plot objs")
    p_list = imap(images, function(img, img_idx){
        img_p_list = imap(features, function(feature, feature_idx){
            p =SpatialPlot(obj, image = img, 
                                features= feature, 
                                pt.size.factor = pt.size.factors[[img_idx]],
                                ...) +
            scale_fill_gradientn(name = feature, 
                                 colours = col_palette,
                                 limits = c(exp_lows[[feature_idx]],
                                            exp_highs[[feature_idx]])
                                 )
            if(img_idx>1){ # remove legend for all plots after first row
                p +  theme(legend.position = "none")
            }else{
                p
            }
        })
        wrap_elements( grid::textGrob(sample_labels[[img_idx]]) ) | img_p_list
    })

    # Decide output format
    if(assemble_col){ # Put all gene rows togehter
        wrap_plots(p_list, nrow = length(images))
    }else{ # Output as a list
        print("Plot not assembled yet. Output as a list")
        p_list
    }
}



SpatialFeaturePlotPlusFeatureAsRow = function(obj, features, images, sample_labels, pt.size.factors, exp_lows, exp_highs, col_palette, assemble_row=T, ...){
    # Seurat plotting function
    message("Creating plot objs")
    p_list = imap(features, function(feature, feature_idx){
        img_p_list = imap(images, function(img, img_idx){
            p =SpatialPlot(obj, image = img,
                                features= feature, 
                                pt.size.factor = pt.size.factors[[img_idx]],
                                ...) +
            scale_fill_gradientn(name = feature, 
                                 colours = col_palette,
                                 limits = c(exp_lows[[feature_idx]],
                                            exp_highs[[feature_idx]])
                                 )
            if(feature_idx==1){ # remove legend for all plots after first row
                p = p + labs(title = sample_labels[[img_idx]]) + 
                        theme(plot.title = element_text(hjust = 0.5, face='bold', size = 15))
                    }
            if(img_idx==1){ # remove legend for all plots after first row
                p = p +  theme(legend.position = "left",
                                legend.title = element_text(size =15))
            }else{
                p = p +  theme(legend.position = "none")
            }
            p
        })
        # wrap_plots(img_p_list, ncol = length(images))
        img_p_list
    })

    # Decide output format
    if(assemble_row){ # Put all gene rows togehter
        wrap_plots(p_list, nrow = length(features))
    }else{ # Output as a list
        print("Plot not assembled yet. Output as a list")
        p_list
    }
    

}

