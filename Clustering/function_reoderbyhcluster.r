PivotWiderToMatrix_ = function(data_use, ident_column, groupby_column, value_column){
        value_mtx = data_use %>% 
        dplyr::select(all_of(c(ident_column,groupby_column, value_column))) %>%
        pivot_wider(id_cols = all_of(ident_column), names_from = all_of(groupby_column), values_from = all_of(value_column), values_fill=0) %>% 
        column_to_rownames(ident_column) %>% 
        as.matrix() %>% 
        t()
    return
}

ReorderByHCluster = function(data_use, ident_column, groupby_column, value_column){
    # 1. convert to wide
    value_mtx = PivotWiderToMatrix_(data_use, ident_column, groupby_column, value_column)
    # 2. cluster using hclust
    hcluster = hclust(dist(value_mtx))
    
    # 3. get the order
    label_ordered = hcluster$labels[hcluster$order]
    
    # 4. reorder the data
    data_use = data_use %>% mutate( {{ groupby_column}} := factor(.data[[groupby_column]], levels = label_ordered))
    message("Reordered groupby_column: ", groupby_column, " by hclust")
    # 5. return the data
    return(data_use)
}


# NOT TESTING YET!  NEED MORE WORK!
ReorderSplitByHCluster = function(data_use, ident_column, groupby_column, value_column, split_by_vector = NULL){
    # No Split
    if(is.null(split_by_vector)) return(ReorderByHCluster(data_use, ident_column, groupby_column, value_column))
    value_mtx = PivotWiderToMatrix_(data_use, ident_column, groupby_column, value_column)
    # Check if split vector has the same lenght as the groupby_column item length
    stopifnot(length(split_by_vector) == length(unique(data_use[[groupby_column]])))
    # For each vector value, if the length > 1, split into smaller matrix and run HClust (column of the matrix)
    splitby_ordred_groupby_ident_vector = table(split_by_vector) %>% imap(function(n_element, split_id){
        # Get the matrix or the selected item
        split_ident_position = split_by_vector == split_id # c(1,1,3,3,5,5) == c(5) -> c(F,F,F,F,T,T)
        sub_mtx = value_mtx[,split_ident_position, drop = FALSE]
        
        # If the length is 1, no need to split, return the ident name as is
        if(ncol(sub_mtx) == 1) return(colnames(sub_mtx))

        # Otherwise hclust
        # 2. cluster using hclust
        hcluster = hclust(dist(value_mtx))
        
        # 3. get the order
        label_ordered = hcluster$labels[hcluster$order]
        return(label_ordered)
    }) %>% unlist %>% as.vector()
    return(splitby_ordred_groupby_ident_vector)
}
