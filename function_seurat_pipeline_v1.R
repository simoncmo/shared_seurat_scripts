# Updated: 2/2021
library(Seurat)
library(tidyverse)
library(future) # For integration anchor mode

# Initial process
spatial_initial_process = function(obj, pca_npcs=30, resolution=0.5, verbose=T, no_normalization =F, method = 'SCT' ){
  # Normalization/PCA
  message(str_glue("Using {pca_npcs} pcs, resolution = {resolution}"))

  if(method == 'SCT'){
    if(no_normalization){
      message('Warning, skipping Normalization step. Make sure object has been normalizated. Might only work for single slice')
      DefaultAssay(obj) ='SCT'
      obj = obj %>% 
        RunPCA(assay = "SCT", verbose = verbose, npcs=pca_npcs) 
    }else{
      message('Using SCTransform for normalization')
      obj = obj %>% 
        SCTransform(assay = "Spatial", return.only.var.genes = FALSE, verbose = verbose) %>% 
        RunPCA(assay = "SCT", verbose = verbose, npcs=pca_npcs)
    }
  }else{ # Using log normal
    if(no_normalization){
      message('Warning, skipping Normalization step. Make sure object has been normalizated')
      DefaultAssay(obj) ='Spatial'
      obj = obj %>% 
        FindVariableFeatures(verbose = verbose) %>% 
        ScaleData(verbose = verbose) %>% 
        RunPCA(verbose = verbose, npcs=pca_npcs) 
    }else{
      message('Using LogNormalize for normalization')
      DefaultAssay(obj) = 'Spatial'
      obj = obj %>%  
        NormalizeData(assay = "Spatial", verbose = verbose) %>% 
        FindVariableFeatures() %>% 
        ScaleData() %>% 
        RunPCA( verbose = verbose, npcs=pca_npcs) 
    }
  }

 # find cluster
    obj = obj %>% FindNeighbors(reduction = "pca", dims = 1:pca_npcs) %>%
        FindClusters(verbose = verbose, resolution = resolution) %>%
        RunUMAP(reduction = "pca", dims = 1:pca_npcs)
    obj
}

# Integration with "merge" or "anchor" mode
spatial_integration = function(st_objs, integration_mode='merge', run_cluster=T, future_worker=50, future_max_gb=40, 
                                RPCA_reference, shared_variable_feature_method='SelectIntegrationFeatures', 
                                cluster_resolution =0.5){
  if(integration_mode=='merge'){
    # Merge integration. Standard integration ----------------------
    print("Using Merge mode")
    ptm = Sys.time()
    st_int = merge(st_objs[[1]], st_objs[-1], add.cell.ids = names(st_objs))
    Sys.time() - ptm; ptm = Sys.time()

    DefaultAssay(st_int) <- "SCT"

    # Variable Features
    print('Using intersection of variable features in all objects')
    shared_variable_feats  = reduce(map(st_objs, VariableFeatures), intersect)
    VariableFeatures(st_int) = shared_variable_feats


  }else if(integration_mode=='anchor'){
    print("Using Anchor mode for batch effect correction. For large object recommand 'RPCA' instead")
    # Anchor based Integration. For sample with batch effect ------------------------
    plan("multiprocess", workers = future_worker) 
    options(future.globals.maxSize = future_max_gb * 1000 * 1024^2) # Increase this if showing "This exceeds the maximum allowed size"
    set.seed(1599) # to replicate
    # select features that are repeatedly variable across datasets for integration
    int_features = SelectIntegrationFeatures(object.list = st_objs)
    st_objs      = PrepSCTIntegration(object.list = st_objs, anchor.features = int_features, verbose = FALSE)
    int_anchor   = FindIntegrationAnchors(object.list = st_objs, normalization.method = "SCT", 
                                          anchor.features = int_features, verbose = FALSE)
    st_int       = IntegrateData(anchorset = int_anchor, normalization.method = "SCT", verbose = FALSE)
    DefaultAssay(st_int) <- "integrated"


  }else if(integration_mode=='RPCA'){
    if(missing(RPCA_reference)){
      stop('For RPCA mode, please specify the index (e.g. c(1,3,8)) of reference obj in the object list using RPCA_reference argument')
    }

    # Using SCT assay instead of NormalizeData/FindVariableFeatures/ScaleData
    int_features <- SelectIntegrationFeatures(object.list = st_objs)
    st_objs      <- map(st_objs, RunPCA, features = int_features)
    int_anchors  <- FindIntegrationAnchors(object.list = st_objs, reference = RPCA_reference, reduction = "rpca", 
        dims = 1:30) 


     # this command creates an 'integrated' data assay
    ptm = Sys.time(); ptm
    st_int       <- IntegrateData(anchorset = int_anchors, dims = 1:30)
    Sys.time() - ptm; ptm = Sys.time()

    # specify that we will perform downstream analysis on the corrected data note that the original
    DefaultAssay(st_int) <- "integrated"

    # Run the standard workflow for visualization and clustering
    st_int <- ScaleData(st_int, verbose = FALSE)

  }else{
    stop("integration_mode only support 'merge', 'anchor' or 'RPCA' mode at the moment")
  }
  
  
  # PCA
  # Assay for RunPCR needs to be "integrated" if using Anchor method. "SCT" for Merge mode
  if(integration_mode=='merge') {
    st_int = RunPCA(st_int, assay = "SCT", verbose = FALSE)
  }else if(integration_mode=='anchor'){
    st_int = RunPCA(st_int, assay = "integrated", verbose = FALSE) 
  }else if(integration_mode=='RPCA'){
    st_int = RunPCA(st_int,, verbose = FALSE) 
  }

  # Cluster
  if(run_cluster){
    print('Clustering the integrated obj')
    st_int %>% 
      FindNeighbors(reduction = "pca", dims = 1:30) %>%
      FindClusters(verbose = FALSE, resolution = cluster_resolution) %>%
      RunUMAP(reduction = "pca", dims = 1:30)
  }else{
    print('Only ran PCA. Still need FindNeighbors/FindClusters/RunUMAP')
    st_int
  }
  
}


# wrapper to run seurat cell type transfer from scNRA/ scRNA to ST 
seurat_cell_type_transfer = function(st_obj, reference_obj, ref_cell_type_col_name='cell_type', 
    output_assay_name = 'cell_type_transfer', assign_dominant=T, output_dominant_col){
   # Check
   if(missing(st_obj)|missing(reference_obj)){stop("Please provide both reference (scRNA/snRNA) and ST object")} 
   if(!'pca' %in% names(st_obj@reductions)){stop("PCA not found in reductions slot of ST obj. Run RunPCA first?")}
   if(!'SCT' %in% names(st_obj@assays)){stop("SCT not found in assay slot of ST obj. Run SCTransform first?")}
   if(!ref_cell_type_col_name %in% names(reference_obj@meta.data)){stop(str_glue("{ref_cell_type_col_name} not found in meta.data slot of reference obj. Check the column name for the cell type col again?"))}
   # Start
   ptm = Sys.time() # timer
   shared_variable_feats = intersect(VariableFeatures(reference_obj), VariableFeatures(st_obj))
    message(glue::glue(
        "Using {length(shared_variable_feats)} shared variable features ",
        "between scRNA and ST to transfer label..."
    ))
    
    Sys.time() - ptm; ptm = Sys.time() # timer
    anchors = FindTransferAnchors(
        reference = reference_obj,
        query = st_obj,
        normalization.method = "SCT",
        features = shared_variable_feats,
        reference.assay = 'SCT',
        query.assay = 'SCT'
    )

    Sys.time() - ptm; ptm = Sys.time() # timer
    # Predict the label transfer
    predictions.assay = TransferData(
        anchorset = anchors,
        refdata = reference_obj@meta.data[, ref_cell_type_col_name],
        prediction.assay = TRUE,
        weight.reduction = st_obj[["pca"]],
        dims = 1:30
    )

    Sys.time() - ptm; ptm = Sys.time() # timer
    print(str_glue("Prediction done. Save result to {output_assay_name} assay"))
    st_obj[[output_assay_name]] = predictions.assay

    # Assign cell type
    if(assign_dominant){
        if(missing(output_dominant_col)){output_dominant_col = str_c(output_assay_name, "_dominant")}
        print(str_glue("Using cell type with max proportions to assign cell type. Save to {output_dominant_col}"))
        st_obj@meta.data[[output_dominant_col]] = getDominantAssignment(st_obj, assay_name=output_assay_name)
    }
    return(st_obj)
}


# Get 
getDominantAssignment=function(obj, assay_name){
  if(!assay_name%in%names(obj@assays)) {stop(str_glue("Assay {assay_name} not found."))}
  celltype_percent = obj@assays[[assay_name]] %>% GetAssayData %>% as.matrix %>% t %>% as.data.frame
  celltype_df     = data.frame(row.names = rownames(celltype_percent),
                               celltype  = colnames(celltype_percent)[max.col(celltype_percent, ties.method = 'first')])
  celltype_df
}
