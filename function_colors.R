# Extrend a color palette with DiscretePalette from Seurat
ExtendColorVector = function(color_vector, target_class, reduce = T, palettes_use = 'glasbey'){
    missing_classes = setdiff(unique(target_class), names(color_vector))
    # Get color and assign
    new_col = Seurat::DiscretePalette(length(missing_classes), palette = palettes_use) %>% setNames(missing_classes)
    # Combine
    color_vector = c(color_vector, new_col)
    # Remove 'Non exisiting item'
    if(reduce) color_vector = color_vector[unique(target_class)]
    return(color_vector)
}
