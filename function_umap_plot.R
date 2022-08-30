## This file define easier plotting of UMAP result
## similar to MC3 package but use ggplot
## need umap
## 8/10/2022

library(umap)
library(tidyverse)

#############################################
# create umap obj
#############################################
MakeUmapObj = function(mtx, meta_data, n_neighbors = 15, seed = 42){
    # Run umap. Note mtx is sample (row) x feature (col)
    set.seed(seed)
    umap_result = umap(mtx, n_neighbors = n_neighbors)
    
    # process meta.data
    if(missing(meta_data)){
        meta_data = data.frame(sample = rownames(mtx), 
                               row.names = rownames(mtx))
    }
    
    # Create obj
    umap_obj = list(
        # Orignal mtx
        data = mtx, 
        # Obj
        coord = umap_result$layout %>% 
            as.data.frame %>% 
            setNames(c('x','y')),
        # Meta.data (meta is sample x metadata(col))
        meta = meta_data[rownames(mtx), ,drop=F]
    )
}

#############################################
# Dim Plot
# Need to be meta in meta table
#############################################

UmapDimPlot = function(umap_obj, group.by, label = F, label_only){
    # Generate data
    p_table = bind_cols(umap_obj$coord, 
                  umap_obj$meta) 
    # Creat plot 
    p = if(!missing(group.by)){
        p_table %>% ggplot(aes(x=x, y=y, color = .data[[group.by]])) 
    }else{
        p_table %>% ggplot(aes(x=x, y=y))
    }
    p = p + geom_point() + 
        cowplot::theme_cowplot() + coord_fixed()
    
    # label 
    if(label){
        # Set what to label
        if(missing(label_only)) label_only = p_table[[group.by]] # if missin label all
        
        # generate label colum
        p_table = p_table %>% mutate(label = ifelse(.data[[group.by]] %in% label_only,
                                                    .data[[group.by]],
                                                    ''
                                                   ))
        
        
        p = p + ggrepel::geom_text_repel(data = p_table, aes(label = label), 
                                        max.overlaps = Inf)
    }
    p
}


############################################################
# Feature PLOT
# DEMO
# gene_plt = c('LGR5','DACH1')
# map(gene_plt, function(gene){
#     # Use External ALL CNV table
#     UmapFeaturePlotCNV(umap_tmp, external_feature = avg_cnv_mtx_all[[gene]], external_feature_name= gene)
# }) %>% wrap_plots()
############################################################
# NOTE EXTERNAL feature allow supplying value NOT in the matrix when doing umap
UmapFeaturePlot = function(umap_obj, feature=NA, palette = RColorBrewer::brewer.pal(9, 'YlOrRd')[3:9], external_feature = NULL, external_feature_name = 'External'){
    # Generate data
    p_table = reduce(umap_obj, bind_cols)
    print(head(external_feature))
    # Creat plot 
    if(!feature %in% colnames(p_table)){
        if(is.null(external_feature)){
            stop(str_glue('Feature {feature} not found in this umap obj. Please or provide external feature'))   
        }else{
            feature = external_feature_name
            p_table[[feature]] = external_feature
        }
    }
    p_table %>% ggplot(aes(x=x, y=y, color = .data[[feature]])) + 
        geom_point() + 
        cowplot::theme_cowplot() + coord_fixed() + 
        scale_color_gradientn(colors = palette) + 
        labs(title = feature)
}

## Special version for CNV since it uses a different color scheme
# NOTE EXTERNAL feature allow supplying value NOT in the matrix when doing umap
UmapFeaturePlotCNV = function(umap_obj, feature=NA, external_feature = NULL, external_feature_name = 'External',
                              palette = c(`0` = "#023858", `0.5` = "#3690C0", `1` = "gray85", `1.5` = "#FEB24C", 
`2` = "#FC4E2A", `2.5` = "#E31A1C", `3` = "#800026")){
    # Process Palette
    cnv_range = if(is.null(external_feature)){
        range(umap_obj$data[[feature]], na.rm = T)
    }else{ # external feature
        range(external_feature, na.rm = T)
    }

    palette_use = palette[seq(cnv_range[[1]], 
                              cnv_range[[2]], by = 0.5) %>% as.character]
    #Plot
    UmapFeaturePlot(umap_obj, feature, palette_use, external_feature, external_feature_name)
}