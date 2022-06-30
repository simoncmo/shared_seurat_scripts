# Extrend a color palette with DiscretePalette from Seurat
ExtendColorVector = function(color_vector, target_class, reduce = T, palettes_use = 'glasbey', palette_method = c('Seurat','RColorBrewer'), randomize = F){
    missing_classes = setdiff(unique(target_class), names(color_vector))
    # Get color and assign
    palette_method = match.arg(palette_method)
    new_col_all = if(palette_method == 'Seurat') Seurat::DiscretePalette(32, palette = palettes_use) else GetRColorBrewerQual()
    
    # rm duplicated color
    new_col_use = setdiff(new_col_all, color_vector) # remove exisiting color
    # Choose color
    new_col_use = if(randomize) sample(new_col_use,length(missing_classes)) else new_col_use[1:length(missing_classes)]
    new_col = new_col_use %>% setNames(missing_classes)
    
    # Combine
    color_vector = c(color_vector, new_col)
    # Remove 'Non exisiting item'
    if(reduce) color_vector = color_vector[unique(target_class)]
    return(color_vector)
}

# Get All Discrete color from RColorBrewer
GetRColorBrewerQual = function(as_list = F){
    df = RColorBrewer::brewer.pal.info %>% 
        filter(category == 'qual') %>% 
        rownames_to_column('palette') 
    palettes = df %>% pmap(function(maxcolors, palette, ...){
     RColorBrewer::brewer.pal(n = maxcolors, name = palette)
     }) %>% setNames(df[['palette']])
    if(as_list) return(palettes) else return(unlist(palettes))
}