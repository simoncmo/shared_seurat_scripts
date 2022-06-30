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