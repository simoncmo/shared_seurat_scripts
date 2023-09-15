
ReorderByHCluster = function(data_use, ident_column, groupby_column, value_column){
    # 1. convert to wide
    value_mtx = data_use %>% 
        dplyr::select(all_of(c(ident_column,groupby_column, value_column))) %>%
        pivot_wider(id_cols = all_of(ident_column), names_from = all_of(groupby_column), values_from = all_of(value_column), values_fill=0) %>% 
        column_to_rownames(ident_column) %>% 
        as.matrix() %>% 
        t()
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
