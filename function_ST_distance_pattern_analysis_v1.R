# For Distance Extraction
# - Created 01/15/2021

# File this function from : PKD/script/PKD_neighborhood_extraction_DWF_clustering_01152021.ipynb

# FUNCTIONS
# General  ----------------------------------------------------------------------------------
# FUNCTIONS
# Get coord_meta df
get_coord_meta_df = function(st, slice, cell_ids, column_include){
    if(missing(slice))    slice    = Images(st)[[1]]
    if(missing(cell_ids)) cell_ids = rownames(st@images[[slice]]@coordinates)
    if(missing(column_include)) column_include = colnames(st@meta.data)
    cbind(st@meta.data[cell_ids, column_include], st@images[[slice]]@coordinates[cell_ids, ])
}

# FUNCTIONS
#functions For distance mtx
get_distance_mtx = function(coord_df, normalize = T){
    dist_mtx = coord_df[,c('imagerow','imagecol')] %>% as.matrix %>% dist %>% as.matrix
    if(normalize){
        min_dist = min(dist_mtx[1,-1])
        dist_mtx = round(dist_mtx / min_dist)
    }
    dist_mtx
}

# List of all neighbor of each spot
get_neighbor_id_list = function(distance_mtx, size = 1){
    # Go through row to get neighbor 
    neighbor_df = apply(distance_mtx, 1, function(distances){
            names(distances[ distances <= size])
    }) %>% setNames(rownames(distance_mtx)) # rowname is the 'center' spot
    
    # Remove self, center spot, from neighbor
    neighbor_list = imap(neighbor_df, ~setdiff(.x , .y))
    neighbor_list
}

# Get neighbor of arbitrary points
get_neighbor_id = function(quary_id, neighbor_list){
    selected_list = neighbor_list[quary_id]
    neighbor_spots = unlist(selected_list) %>% setdiff(quary_id)
    neighbor_spots
}

# Function
# From target spot, get source spot and distance 
get_source_spots = function(source_spots, neighbor_spots, distance_mtx, include_source_spots=T){
    # Id of source
    source_position    = distance_mtx[source_spots, neighbor_spots] %>% apply(2, which.min) # Value is Source's row in table
    # Distance to source 
    distance_vector    = distance_mtx[source_spots, neighbor_spots] %>% apply(2, min)
    # Combine
    source_dist_df     = data.frame(neighbor_spots = names(source_position), # Name is target spots
                                    source_spot = spots_of_interest[source_position],
                                    distance_to_source    = distance_vector[names(source_position)])  # Keep in same order
    # Add source spots
    if(include_source_spots){
        source_to_source_df = data.frame(neighbor_spots = source_spots, 
                                         source_spot = source_spots,
                                         distance_to_source = 0)
        rownames(source_to_source_df) = source_spots
        source_dist_df = rbind(source_dist_df, source_to_source_df)
    }
    source_dist_df
}

# Add meta of source 
add_source_meta = function(meta_df, source_dist_df, meta_columns = c('histology'), add_column_name = 'source'){
    # Extract source meta info
    # Assuming rownames is the 'spot_id'
    source_meta_df = meta_df[,meta_columns, drop=F] %>% 
        setNames(str_c(add_column_name, names(.), sep='_')) %>% 
        rownames_to_column('source_spot')
    
    # Add to source_dist_df 
    source_dist_df = left_join(source_dist_df, source_meta_df, by = 'source_spot')
    return(source_dist_df)
}

# function
# Wrapper 
# # 3. Get final annotated neighborhood spot df
# # With Coordinates, source, distance, source meta, neighborhood spot meta
get_neighborhood_df = function(source_spots, neighboor_spots, distance_mtx, include_source_spots=T,
                    coord_meta_df, 
                    source_meta_columns = c('histology'), 
                    #neighboor_meta_columns = c('histology'),
                    add_source_name = 'source',
                    add_neighbor_name = 'neighbor'){
    # Get 'source' spot and their meta for each gene of interests 
    source_dist_df = get_source_spots(source_spots, neighboor_spots, distance_mtx, include_source_spots) 
    source_dist_df = add_source_meta(coord_meta_df, source_dist_df, meta_columns = source_meta_columns, 
                                     add_column_name = add_source_name)
    
    # Create final annotation df
    coord_meta_df = coord_meta_df %>% 
                setNames(str_c(add_neighbor_name, names(.), sep='_')) %>% 
                rownames_to_column('neighbor_spots') # Target (neighbor) point meta 
    neighborhood_annotated_df = right_join(coord_meta_df, 
                                           source_dist_df, # source point meta and distance 
                                           by = 'neighbor_spots')
    neighborhood_annotated_df
}

# function
# Add features
add_neighboor_features = function(st_obj, neighborhood_df, features, assay){
    if(!missing(assay)) DefaultAssay(st_obj) = assay
    # Spots 
    neighbor_spots = neighborhood_df[, 'neighbor_spots']
    feature_df = FetchData(object = st_obj, vars = features, cells = neighbor_spots)[neighbor_spots,,drop=F] # Get feature for st object
    return(cbind(neighborhood_df, feature_df)) 
}


# DE function  ----------------------------------------------------------------------------------

# DE function
FindNeighborhoodDE = function(obj, neighborhood_df, source_range = 0,
                              distant_spot_range = 3
                             ){

    source_range_spots  = neighborhood_df %>% filter(distance_to_source %in% source_range) %>% pull(neighbor_spots)
    distant_range_spots = neighborhood_df %>% filter(distance_to_source %in% distant_spot_range) %>% pull(neighbor_spots)

    FindMarkers(obj = obj[['SCT']], cells.1 = source_range_spots, 
               cells.2 = distant_range_spots, ) %>% rownames_to_column('gene')
}

SweepNeighborhoodDE = function(obj, neighborhood_df, source_range=0, sweep_range=3:5, sweep_size=1, sweep_overlap =T){
    # Remove range overlap with source 
    sweep_range = setdiff(sweep_range, source_range) 
    
    # Set sweep range
    if(sweep_overlap){
    # Overlap
        sweep_list = map(sweep_range, ~seq(from = ., length.out = sweep_size))
    }
    else{
        # No Overlap
        sweep_group_min = seq(from = min(sweep_range), to = max(sweep_range), by = sweep_size)
        sweep_list = map(sweep_group_min , ~seq(from = ., length.out = sweep_size))
    }

    # Run list of range
    map(sweep_list, function(sweep_group){
        FindNeighborhoodDE(obj, neighborhood_df, source_range, distant_spot_range = sweep_group) %>% 
        mutate(source_range = toString(source_range),
               distant_range = toString(sweep_group)
              )
    }) %>% bind_rows()
}

GetSweepNeighborhoodDEJaccardScore = function(neighborhood_de_sweep_df, fdr_cutoff = 0.05, complete_table =T){
    # Get gene list split by sweep groups
    sweep_gene_list = neighborhood_de_sweep_df %>% 
                            filter(p_val_adj < fdr_cutoff) %>% 
                            group_by(distant_range) %>% 
                            split(f = .$distant_range, x=.$gene)

    # Get Jaccard score
    df = combn(names(sweep_gene_list), m = 2) %>% t %>% as.data.frame %>% 
    pmap(function(V1, V2){
        intersect_gene = intersect(sweep_gene_list[[V1]], sweep_gene_list[[V2]])
        union_gene = union(sweep_gene_list[[V1]], sweep_gene_list[[V2]])
        jaccard = length(intersect_gene)/length(union_gene)
        data.frame(Group1 =V1, 
                   Group2 =V2,
                   Group1_n = length(sweep_gene_list[[V1]]), 
                   Group2_n = length(sweep_gene_list[[V2]]),
                   intersect_n = length(intersect_gene), 
                   union_n = length(union_gene), 
                  Jaccard_score = jaccard)
    }) %>% bind_rows()
    
    if(complete_table){
        df = rbind(df, 
                   df %>% dplyr::rename(Group2 = Group1, Group1= Group2)  %>% 
                        dplyr::rename(Group2_n = Group1_n, Group1_n= Group2_n)
                  )
    }
    df
}

# PLOT FUNCtion
## Expression in the neighborhood
plot_spatial_neighboor_feature = function(neighborhood_df, 
                                          x_col = 'neighbor_imagecol',
                                          y_col = 'neighbor_imagerow',
                                          flip_axis = c('none','x','y','both'),
                                          pt_size = 1.5,
                                          feature,
                                          rev_color = F,
                                          palette
                                         ){
    if(missing(feature)) stop('Missing feature to plot')
    # Adjust axis
    flip_axis = match.arg(flip_axis)
    if(flip_axis %in% c('x','both')) neighborhood_df[[x_col]] = -neighborhood_df[[x_col]]
    if(flip_axis %in% c('y','both')) neighborhood_df[[y_col]] = -neighborhood_df[[y_col]]
    
    # Color
    col_palette = Seurat:::SpatialColors(10)
    if(rev_color) col_palette = rev(col_palette)
    if(!missing(palette)) col_palette = palette

   neighborhood_df %>% 
    ggplot(aes(x = .data[[x_col]], y = .data[[y_col]], color = .data[[feature]]))+ 
        geom_point(size = pt_size) + 
        coord_fixed() + 
        scale_color_gradientn(colors = col_palette) +
        labs(title = feature)+
        theme_void() +
        theme(plot.title = element_text(size = 20, face = 'bold', hjust = 0.5))
}

# Function
# Distance. A wrapper
plot_spatial_neighboor_distance = function(neighborhood_gene_df, 
                                           distance_col_name = 'distance_to_source',
                                           plot_title = 'Distance to Source',
                                           ...){
    loc_palette = rev(viridisLite::inferno(10))#paletteer::paletteer_d("colorBlindness::Blue2DarkOrange12Steps")#paletteer::paletteer_d("awtools::a_palette") #RColorBrewer::brewer.pal(n=9, 'BuPu')
    plot_spatial_neighboor_feature(neighborhood_gene_df, feature=distance_col_name, palette = loc_palette, ...) + 
    labs(title = plot_title)
}

# Wrapper
# Plot multiple features and add distance 
plot_spatial_neighboor_features = function(neighborhood_gene_df, 
                                           features, 
                                           show_distance_plot = T,
                                           distance_col_name = 'distance_to_source', ...){
    # Feature plots
    plist = map(features, function(feature){
       plot_spatial_neighboor_feature(neighborhood_gene_df, feature = feature, ...) 
    }) %>% setNames(features)
                
    # Add distant plot
    if(show_distance_plot){
        distance_p = plot_spatial_neighboor_distance(neighborhood_gene_df, distance_col_name, ...)
        plist = c(list(distance = distance_p), plist)
        }
    # Combine plot
    wrap_plots(plist)
}

plot_feature_distance_curve = function(neighborhood_feature_df,
                                 features,
                                 group.by = 'source_histology',
                                 idents,
                                 distance_col = 'distance_to_source',
                                 data_type = 'Gene',
                                 value_type = 'Expression', # use 'Percentage or probability' for Cell type
                                 smooth = T
                                ){
    if(missing(idents)) idents = unique(neighborhood_feature_df[[group.by]])
    # Facet formula
    facet_formula = str_c(data_type,'~', group.by) %>% as.formula

    cols_to_keep = c(group.by, distance_col, features)
    df = neighborhood_feature_df %>% 
        filter(.data[[group.by]] %in% idents) %>% 
        select(all_of(cols_to_keep)) %>% a
        pivot_longer(cols = -c(.data[[distance_col]], .data[[group.by]]), 
                     names_to = data_type, 
                     values_to = value_type) 
    df_Avgexp = df %>% group_by(.data[[distance_col]], .data[[group.by]], .data[[data_type]]) %>% 
        summarise(Ave_exp = mean(.data[[value_type]], na.rm=T))
    
    if(smooth){ # Use smooth line
        p = ggplot(df, aes(x = .data[[distance_col]], 
               y = .data[[value_type]], 
               color = .data[[group.by]]) ) + 
        facet_grid(facet_formula)+
        geom_smooth(aes(fill = .data[[group.by]]), method = 'loess', formula = 'y~x', alpha = 0.15, size = 0, span=0.75) + # SE. span control smoothness
        stat_smooth(geom="line", method = 'loess', formula = 'y~x', alpha=0.6, size=2, span = 0.75) #span control smoothness
    } else { # Use Ave expression line plot
        p = ggplot(df_Avgexp, aes(x = .data[[distance_col]], 
               y = Ave_exp, 
               color = .data[[group.by]]) ) + 
        facet_grid(facet_formula)+
        geom_line() + 
        geom_point()
    }
    p + theme_bw() + 
        theme(panel.grid = element_blank(),
              strip.background = element_rect(fill = 'gray90', color = 'transparent')
                 )
}

