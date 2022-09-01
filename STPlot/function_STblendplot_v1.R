library('Seurat')

#source('/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/NMK/NMK_integration_v5/script/function_ST_plot_general_v1.R')


# Main function
STBlendPlot = function(st, feature_pairs, images=1, show_histology=F, show_blend_only = F, wrap_slices=T,...){
  if(typeof(feature_pairs)!='list') stop('STOP: feature_pairs need to be a list')
  #show_blend_only   =ifelse(length(feature_pairs)>1, T, F) 
  MultiFeatureSTBlendPlot(st, feature_pairs, images, show_blend_only=show_blend_only, show_histology=show_histology, wrap_slices=T,...)
}

SThistology = function(st, images){
  if(is.numeric(images)) images = Images(st)[images]
  p = SpatialPlot(st, images = images, pt.size.factor=0, stroke=NA, combine=T) 

  p & scale_x_continuous(expand = c(0,0)) &
      scale_y_continuous(expand = c(0,0)) &
      theme(plot.title = element_text(face = 'bold', hjust=0.5),
            plot.margin = unit(c(0,0,0,0), 'lines')) &
      NoLegend()
}

# Make plot for multiple genes for multiple slice
## Need to set limits for color scalle
MultiFeatureSTBlendPlot = function(st, feature_pairs, images=1, show_blend_only, show_histology=F, wrap_slices=T,...){
  plist = map(feature_pairs, function(feature_pair){
    p_slices = MultiImageSTBlendPlot(st, feature_pair=feature_pair, images=images, combine=F, show_blend_only=show_blend_only, show_histology=show_histology,...)
    if(wrap_slices) p_slices %>% wrap_plots(nrow=1)
    else p_slices
  }) 
  # histology
  # if(show_histology){
  #   histology = SThistology(st, images = images)
  #   return(wrap_plots(c(list(histology), plist), ncol=1))
  # }else{
  #   plist %>% wrap_plots(ncol=1) 
  # } 
  if(!wrap_slices) return(plist)
  plist %>% wrap_plots(ncol=1) 
}

# Make plot for 1 genes for multiple slice
MultiImageSTBlendPlot = function(st, feature_pair, images=1, combine=F, assay = 'SCT', show_blend_only, show_histology=F, 
                                  flip_x_images=NA, flip_y_images=NA, # use to flip images
                                  ...){
  # Select images 
  imgs   = if(is.numeric(images)) Images(st)[images] else images


  # Get cell id for normalizing 
  normalizing_cellids = map(images, function(img){
    rownames(st@images[[img]]@coordinates)
  }) %>% unlist

  # Get exp data from all slices to speed up 
  DefaultAssay(st) = assay
  exp_df_all = FetchData(st, vars = feature_pair, slot = 'data', cells = normalizing_cellids)

  # Plot
  map(imgs, function(img){
    flip_x = ifelse(img %in% flip_x_images, T,F)
    flip_y = ifelse(img %in% flip_y_images, T,F)
    plist = SingleSTBlendPlot(st, features = feature_pair, image = img, combine = combine, assay = assay,
      show_blend_only=show_blend_only, 
      normalizing_cellids = normalizing_cellids,
      multi_slice_exp_df  = exp_df_all,
      flip_x=flip_x, flip_y=flip_y,
     ...) 
    if(show_histology){
      histology = SThistology(st, images = img) + labs(title = img)
      p = wrap_plots(c(list(histology), plist), ncol=1)
    }else{
      p = plist %>% wrap_plots(ncol=1)
    }
    p
  }) %>% wrap_plots(nrow=1)

}


### V3 Current working version 8/3/2021
# Adding new scaling function 8/31/2022
SingleSTBlendPlot = function(st, features, combine = T, assay = 'SCT', image, show_blend_legend =T,
                       normalizing_cellids, normalize_among_slices = T, # This is used when plotting multiple slices
                       multi_slice_exp_df, # Used for multiple slice and speed up fetching exp data
                       show_blend_only= F, plot_title = c('gene','slice','slice-gene'),
                       flip_x =F, flip_y=F,
                       pt.size.factor=1,
                       feature_colors =c("#ff0000", "#00ff00"), negative_color = 'gray10',
                       # Allow x, y scaling : 8/31/2022
                       scale_x = 1, scale_y = 1
                            ){
  # Select cells for normalization
  if(missing(normalizing_cellids) | !normalize_among_slices ) normalizing_cellids = rownames(st@images[[image]]@coordinates)

  # Get data
  DefaultAssay(st) = assay
  exp_df = if(!missing(multi_slice_exp_df)){
                multi_slice_exp_df[normalizing_cellids, ]
              }else{
                FetchData(st, vars = features, slot = 'data', cells = normalizing_cellids)
              }
  
  # Fill NA
  exp_df[is.na(exp_df)] = 0
  
  # Blend: Normalize to 0-9, then just feature1 (0-9) + feature2 (0-9) * 10
  exp_blend = blendExpression_v2(exp_df, features)

  # Keep only cells in the image for plotting
  exp_blend = exp_blend[rownames(st@images[[image]]@coordinates), ]
  
  # Legend
  col_blend = Seurat:::BlendMatrix(negative.color = negative_color, two.colors = feature_colors)
  
  # Location
  if(missing(image)){ image = Images(st)[1] }
  exp_df_blend = cbind(st@images[[image]]@coordinates, exp_blend) %>% mutate(across(.cols=everything(), .fns=as.numeric)) 
  
  # Center the plot
  exp_df_blend = exp_df_blend %>% st_scale_and_center()
  
  # titles
  if(length(plot_title)>1) plot_title = 'slice-gene'
  title_list = list()
  title_list$`1` = switch(plot_title, 'gene' = features[1], 'slice' = image, 'slice-gene' = str_glue('{image} {features[1]}'))
  title_list$`2` = switch(plot_title, 'gene' = features[2], 'slice' = image, 'slice-gene' = str_glue('{image} {features[2]}'))
  title_list$`3` = switch(plot_title, 
                          'gene' = str_c(features[1],features[2], sep='-'), 
                          'slice' = image, 
                          'slice-gene' = str_glue('{image} {str_c(features[1],features[2], sep="-")}'))
        
  # New Aug 2022 - allow scaling of x and y axis
  exp_df_blend = exp_df_blend %>% mutate(imagecol = imagecol * scale_x)
  exp_df_blend = exp_df_blend %>% mutate(imagerow = imagerow * scale_y)
  
  # Plot - changed to imagecol/imagerow 8/31/2022
  plist = list()  
  if(!show_blend_only){
    plist[[features[[1]] ]] = exp_df_blend %>%
      ggplot(aes_string(x = 'imagecol', y='imagerow', color = paste0('`',features[[1]],'`') )) + 
      scale_color_gradientn(colors = col_blend[,1], limits = c(0,9)) + # Set color range be 0-9 to get correct color range
      labs(title = title_list[1])
    
    plist[[features[[2]] ]] = exp_df_blend %>%
      ggplot(aes_string(x = 'imagecol', y='imagerow', color = paste0('`',features[[2]],'`'))) + 
      scale_color_gradientn(colors = col_blend[1,], limits = c(0, 9)) + # Set color range be 0-9 to get correct color range
      labs(title = title_list[2])
  }
  # Blend plot
  plist[[str_c(features, collapse ="_")]] = exp_df_blend %>%
    ggplot(aes_string(x = 'imagecol', y='imagerow', color = paste0('`',str_c(features, collapse ="_"),'`')   )) + 
    scale_color_gradientn(colors = as.vector(col_blend), limits = c(0,99)) + # Set color range be 0-99 to get correct color range
    labs(title = title_list[3])

  # axis orientation - changed 
  #x_axis = if(flip_x) scale_x_reverse(expand = c(0,0), limits = c(64,-64)) else scale_x_continuous(expand = c(0,0), limits = c(-64,64))
  #y_axis = if(flip_y) scale_y_reverse(expand = c(0,0), limits = c(64,-64)) else scale_y_continuous(expand = c(0,0), limits = c(-64,64))
  x_axis = if(flip_x) scale_x_reverse(expand = c(0,0)) else scale_x_continuous(expand = c(0,0))
  y_axis = if(flip_y) scale_y_reverse(expand = c(0,0)) else scale_y_continuous(expand = c(0,0))
 
  
  # Apply shared setups
  plist = map(plist, ~. + geom_point(size=pt.size.factor) + coord_fixed() + theme_void() + 
                x_axis +
                y_axis + 
                theme(plot.title = element_text(face = 'bold', hjust=0.5),
                      plot.margin = unit(c(0,0,0,0), 'lines')
                     ) + 
                NoLegend() )

  # Set to 1 color if all 0
  plist = imap(plist, function(p, idx){ 
    if(max(exp_blend[,idx]) == 0){
      suppressMessages(p + scale_color_gradientn(colors = col_blend[1,1])) 
    }else{
      p
    }
  })
  
  # Legend
  if(show_blend_legend){
    plist$legend = BlendMap_v2(col_blend, exp_df_blend) 
    plist$legend = suppressMessages(plist$legend + 
                                      scale_x_continuous(breaks = quantile(0:9, seq(0.1, 1, length.out=4)), 
                                                         labels = quantile(exp_df[,1], seq(0.1 ,1,length.out=4)) %>% round(1),
                                                         expand = c(0,0)) + 
                                      scale_y_continuous(breaks = quantile(0:9, seq(0.1 , 1, length.out=4)), 
                                                         labels = quantile(exp_df[,2], seq(0.1 , 1,length.out=4)) %>% round(1),
                                                         expand = c(0,0))
    ) + coord_fixed() + theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), 'lines'))
  }
  
  # Combine
  if(combine){
    return(wrap_plots(plist))
  }
  return(plist)
}




## Normalize
normalize_01 = function(x, to_round=F, to_int=F){ # Normalize a vector of number to 0 to 1 range
  x = if(all(is.na(x))| all(x %in% 0)){ rep(0, length(x)) } # If All NA or 0. Set all to 0
      else{
        (x-min(x))/(max(x)-min(x))
      }
  if(to_round) x = round(x, 1)
  if(to_int)   x = round(9*x)
  x
}

## Blend exp
blendExpression_v2 = function(exp_df, features){
  exp_blend = exp_df %>% mutate(across(everything(), normalize_01, to_round = T, to_int = T))
  exp_blend[,3] = exp_blend[,1] + exp_blend[,2] * 10 
  colnames(exp_blend) <- c(features, paste(features, collapse = "_"))
  return(exp_blend)
}



## Legend
BlendMap_v2 = function (color.matrix, coord_exp_df) 
{
  coord_names = c("tissue", "row", "col", "imagerow", "imagecol")
  exp_df = coord_exp_df[,!(names(coord_exp_df) %in% coord_names), drop=F] # Drop =F keep data.frame even only have 1 column
  row_name = colnames(exp_df)[1]
  col_name = colnames(exp_df)[2]
  xrange = range(as.numeric(exp_df[,1]))
  yrange = range(as.numeric(exp_df[,2]))
  
  color.heat <- matrix(data = 1:prod(dim(color.matrix)) - 1, 
                       nrow = nrow(color.matrix), 
                       ncol = ncol(color.matrix), 
                       dimnames = list(1:nrow(color.matrix), 
                                       1:ncol(color.matrix)))
  xbreaks <- seq.int(from = 0, to = nrow(color.matrix), by = 2)
  ybreaks <- seq.int(from = 0, to = ncol(color.matrix), by = 2)
  xlabels =  seq(from = xrange[1], to = xrange[2], length.out = length(xbreaks))
  ylabels =  seq(from = yrange[1], to = yrange[2], length.out = length(ybreaks))
  color.heat <- Seurat:::Melt(color.heat)
  color.heat$rows <- as.numeric(as.character(color.heat$rows))
  color.heat$cols <- as.numeric(as.character(color.heat$cols))
  color.heat$vals <- factor(color.heat$vals)
  plot <- ggplot(data = color.heat, mapping = aes_string("rows", 
                                                         "cols", fill = "vals")) + 
    geom_raster(show.legend = FALSE) + 
    theme(plot.margin = unit(rep.int(0, times = 4), units = "cm")) + 
    scale_x_continuous(breaks = xbreaks, expand = c(0, 0), labels = xlabels) + 
    scale_y_continuous(breaks = ybreaks, expand = c(0, 0), labels = ylabels) + 
    scale_fill_manual(values = as.vector(color.matrix)) + 
    labs(x = row_name, y = col_name) +
    cowplot::theme_cowplot()
  return(plot)
}




# Use adjust data in the FeaturePlot
adjust_data = function(data, cols){
  features <- colnames(x = data)[4:ncol(x = data)]
  cells    = rownames(data)
  min.cutoff <- map(features, ~min(data[, .]))
  max.cutoff <- map(features, ~max(data[, .]))
  check.lengths <- unique(x = vapply(X = list(features, min.cutoff, 
                                              max.cutoff), FUN = length, FUN.VALUE = numeric(length = 1)))
  brewer.gran <- ifelse(test = length(x = cols) == 1, yes = RColorBrewer::brewer.pal.info[cols, 
  ]$maxcolors, no = length(x = cols))
  
  if (length(x = check.lengths) != 1) {
    stop("There must be the same number of minimum and maximum cuttoffs as there are features")
  }
  
  scale_fun =  function(index) {
    data.feature <- as.vector(x = data[, index])
    min.use <- Seurat:::SetQuantile(cutoff = min.cutoff[index - 3], data.feature)
    max.use <- Seurat:::SetQuantile(cutoff = max.cutoff[index - 3], data.feature)
    data.feature[data.feature < min.use] <- min.use
    data.feature[data.feature > max.use] <- max.use
    
    data.cut <- if (all(data.feature == 0)) {
      0
    }
    else {
      as.numeric(x = as.factor(x = cut(x = as.numeric(x = data.feature), breaks = brewer.gran)))
    }
    return(data.cut)
  }
  
  data[, 4:ncol(x = data)] <- sapply(X = 4:ncol(x = data), scale_fun)
  
  
  colnames(x = data)[4:ncol(x = data)] <- features
  rownames(x = data) <- cells
  return(data)
}
