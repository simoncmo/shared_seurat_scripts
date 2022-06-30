SplitSeuratObj = function(obj, split.by, use_future = F){
    if(missing(split.by)) stop("Need identity to split object")
    if(! split.by %in% colnames(obj@meta.data)) stop(str_glue("{split.by} identity not found in thie object"))
    
    Idents(obj) = split.by
    all_idents  = unique(obj@meta.data[[split.by]])
    
    # Future
    if(use_future){
        message('Uses Future for faster subsetting')
        future::plan("multicore", workers = future::availableCores())
        options(future.globals.maxSize = 100*1024^3, future.rng.onMisuse= 'ignore')
    }
    map_FUN  = if(use_future) furrr::future_map else map

    obj_list = map_FUN(all_idents, function(ident){
        subset(obj, idents = ident)
    }) %>% setNames(all_idents)
}