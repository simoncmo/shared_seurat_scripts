GetNeighborSpots = function(st, target_ids, distance = 4, min_distance = 0, slice, include_target = F){
    if(missing(st)|missing(target_ids)) stop('missing ST obj or Target ids')
    if(missing(slice)) slice = Images(st)[[1]]
    
    # Get distance matrix and normalize
    min_dist = st@images[[slice]]@coordinates %>% select(row, col) %>% dist(diag = F, upper=T) %>% min
    dist_mtx = st@images[[slice]]@coordinates %>% select(row, col) %>% dist(diag = T, upper=F) %>% as.matrix 
    dist_mtx = round(dist_mtx/min_dist)
    
    # Get IDS
    neighbor_ids = dist_mtx[,target_ids] %>% apply(1, min) %>% .[.<=distance] %>% .[.>=min_distance] %>% names
    neighbor_ids = if(include_target) neighbor_ids else setdiff(neighbor_ids, target_ids)
    
    return(neighbor_ids)
}

################
# GET SPATIAL RATIO GIVEN SPOTS
################
# API: 
# Usage : PlotBarplotCellRatioPerArea(obj = STobj, assay_name = 'RCTDcell', area_list = dbscan_area_list)
PlotCellRatioBarplotPerArea = function(obj, assay_name, area_list, cell_include, 
                                       cell_palette, label =T, label_threshold = 0.05, 
                                       label_col = 'black', show_percent = T){
    df = GetCellRatioPerArea(obj, assay_name, area_list, cell_include)
    CellRatioBarplotPerArea(df, cell_palette, label, label_threshold, 
                                       label_col, show_percent)
}

# Usage: PlotBarplotCellRatioByIdent(obj = STobj, assay_name = 'RCTDcell', group.by = 'seurat_clusters')
# This API uses Identity in Seurat meta.data to define cell type to plot
PlotBarplotCellRatioByIdent = function(obj, assay_name, group.by, ident_include, cell_include, 
                                       cell_palette, label =T, label_threshold = 0.05, 
                                       label_col = 'black', show_percent = T){
    # Grab identity 
    area_list= GetSpotFromIdent(obj, group.by = group.by, ident_include=ident_include)
    
    # Get DF
    df = GetCellRatioPerArea(obj, assay_name, area_list, cell_include)
    # PLOT
    CellRatioBarplotPerArea(df, cell_palette, label, label_threshold, 
                                       label_col, show_percent)
}

# ------------ INTERNAL --------------
## COLOR
GetDiscreteColors = function(method = c('RColorBrewer','Seurat'), n){
    method = match.arg(method)
    color_pal = switch(EXPR = method, 
           RColorBrewer = GetAllBrewerColor(),
          Seurat = Seurat::DiscretePalette(n = 36))
    color_pal = if(!missing(n)) color_pal[1:n] else color_pal
    return(color_pal)
}
GetAllBrewerColor = function(palette_type = 'qual'){
    pmap(RColorBrewer::brewer.pal.info %>% filter(category == all_of(palette_type)) %>% rownames_to_column('palette'), function(maxcolors, palette, ...){
    RColorBrewer::brewer.pal(n = maxcolors, name = palette)
    }) %>% reduce(c)
}

## SPOT
GetSpotFromIdent = function(obj, group.by, ident_include){
    # Select IDENT
    idents_to_plot = obj@meta.data[[group.by]] %>% unique
    if(!missing(ident_include)) idents_to_plot = idents_to_plot %in% ident_include
    
    # Make list
    map(idents_to_plot, function(ident){
        obj@meta.data %>% filter(.data[[group.by]] %in% all_of(ident)) %>% rownames
    }) %>% setNames(idents_to_plot)
}

GetCellRatioPerArea = function(obj, assay_name, area_list, cell_include, pivot = T){
    if(any(c(missing(obj), missing(assay_name), missing(area_list)))) stop('missing essential arguments')
    if(missing(cell_include)) cell_include = rownames(obj@assays[[assay_name]]) # use all cells
    
    percent_per_area_list = map(area_list, function(ids){
        obj@assays[[assay_name]] %>% GetAssayData %>% as.matrix %>% 
        .[cell_include, ids, drop=F] %>% # Select spots
        rowSums() %>% # Sum through all spots
        (function(x){x/sum(x)})(.) # Normalize to 1
    })
    
    outdf = percent_per_area_list %>% 
            bind_rows() %>% t %>% as.data.frame %>% 
            setNames(names(area_list)) %>% 
            rownames_to_column('CellType')
    if(pivot) outdf = outdf %>% pivot_longer(cols = -CellType, names_to = 'Area', values_to = "Percentage")
    return(outdf)
}

CellRatioBarplotPerArea = function(df, cell_palette, label =T, label_threshold = 0.05, label_col = 'gray90', show_percent = T){
    celltypes = df$CellType %>% unique
    if(missing(cell_palette)) cell_palette = GetDiscreteColors(method = 'RColorBrewer') 
    
    # Label 
    df = df %>% mutate(Label = ifelse(Percentage > label_threshold, CellType, ''))
    if(show_percent) df = df %>% mutate(Label = ifelse(Percentage > label_threshold, 
                                                       str_glue("{Label}\n({round(Percentage,2)*100}%)"),
                                                       ''))
    # MainPlot
    p = ggplot(df, aes(x = Area, y = Percentage , fill = CellType)) + 
        geom_bar(stat = 'identity') + 
        labs(title = 'Cell Type Percentage in each Area') +
        scale_fill_manual(values = cell_palette) +
        cowplot::theme_cowplot()
    
    # Label
    p = if(!label) p else p + geom_text(aes(label = Label), position = position_stack(vjust = 0.5), color = label_col)
    
    return(p)
}

