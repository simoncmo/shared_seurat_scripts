library(tidyverse)
library(Seurat)

message("Functions includes: FindMarkersFormated, FindMarkersEachVsRef, FindMarkersEachVsRefComplete")
message("FindMarkers with better format: FindMarkersFormated")
FindMarkersFormated = function(obj, group.by = 'seurat_clusters', ident.1=NULL, ident.2=NULL, ...){
    message("group.by use: ", group.by)
    Idents(obj) = group.by
    message("Ident1 names as cluster to be same as FindAllMarkers")
    FindMarkers(obj,ident.1 = ident.1, ident.2=ident.2,  ...) %>% 
        rownames_to_column('gene') %>%
        mutate(cluster = as.character(ident.1),
            ident2 = ident.2,
            group_by = group.by
        ) 
}

# Compare every other ident in the group.by to ref
message("Compare every other ident in the group.by to ref: FindMarkersEachVsRef")
FindMarkersEachVsRef = function(obj, group.by='seurat_clusters', ident_ref, ...){
    message('Ref. ident is: ', ident_ref)
    remain_idents = setdiff(obj@meta.data[[group.by]], ident_ref)
    message('Remain idents are: ', toString(remain_idents))
    map(remain_idents, function(ident_1_use){
        message('Ident1 use: ', ident_1_use)
        FindMarkersFormated(obj, group.by = group.by, ident.1 = ident_1_use, ident.2 = ident_ref, ...)
    }) %>% bind_rows()
}

# Compare every other ident in the group.by to ref
# Also include Ref vs NotRef
message("Compare every other ident in the group.by to ref, also include Ref vs NotRef: FindMarkersEachVsRefComplete")
FindMarkersEachVsRefComplete = function(obj, group.by='seurat_clusters', ident_ref, ...){
    message('Ref. ident is: ', ident_ref)
    # 1. Ref vs Not Ref
    deg_ref_notRef_df = FindMarkersFormated(ST_use, group.by = st_groupby, ident.1 = ident_ref) %>% 
        mutate(ident2 = str_glue("Not_{ident_ref}"))
    
    # 2. Each other ident vs Ref
    remain_idents = setdiff(obj@meta.data[[group.by]], ident_ref)
    message('Remain idents are: ', toString(remain_idents))
    deg_each_vs_ref_df = map(remain_idents, function(ident_1_use){
        message('Ident1 use: ', ident_1_use)
        FindMarkersFormated(obj, group.by = group.by, ident.1 = ident_1_use, ident.2 = ident_ref, ...)
    }) %>% bind_rows()

    # 3. Combine
    return(bind_rows(deg_ref_notRef_df, deg_each_vs_ref_df))
}