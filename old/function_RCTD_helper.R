# Helper function for RCTD package
library(Seurat)
library(tidyverse)
library(RCTD)

equal_downsample = function(obj, n =1000, cell_type_col = "liger_ident_coarse"){
  meta = data.frame(barcode = colnames(obj),
                    cell_type = obj@meta.data[[cell_type_col]] )
  # sample per cell type
  # id_keep = meta %>% group_by(cell_type) %>% sample_n(size = n) %>% .$barcode
  
  id_keep = lapply(unique(meta$cell_type), function(type){
    tmp_df = meta %>% filter(cell_type == type)
    tmp_df$barcode %>% sample(size = min(nrow(tmp_df), n))
  }) %>% unlist()
  
  # subset Seurat obj
  obj = obj[,id_keep]
  obj        
}

# Help for scRNA processing - To get scRNA obj directly from Seurat obj
dgeToSeurat_direct = function (dge, meta_data, cell_type_dict, save_rds=F, refdir='Data/Reference')
{
  # DGE
  raw.data <- as(dge, "dgCMatrix")
  rownames(raw.data) <- rownames(dge)
  
  # Meta
  rownames(meta_data) = meta_data$barcode
  meta_data$barcode = NULL
  
  # Shared barcode filtering
  common_barcodes = intersect(colnames(raw.data), rownames(meta_data))
  raw.data = raw.data[, common_barcodes]
  meta_data = meta_data[common_barcodes, ]
  
  # Cell type info
  cell_type_dict$Name <- factor(cell_type_dict$Name)
  rownames(cell_type_dict) = cell_type_dict[, "Cluster"]
  cell_ident = meta_data$cluster
  true_type_names = lapply(cell_ident, function(x) cell_type_dict[as.character(x), "Name"])
  true_type_names = unlist(true_type_names)
  meta_data$liger_ident_coarse = true_type_names
  
  
  # reference file
  reference = Seurat::CreateSeuratObject(raw.data, meta.data = meta_data)
  
  # Save RDS file 
  if(save_rds){
    dir.create(refdir, recursive = T)
    saveRDS(reference, paste(refdir, "SCRef.RDS", sep = "/"))
    dref <- RCTD:::create_downsampled_data(reference, refdir) # Might not work
  }else{
    # dref <- RCTD:::create_downsampled_data(reference, save.file=F) # Not working not sure why
    dref <- equal_downsample(reference, n = 1000)
  }
  
  return(dref)
}

RCTD_ref_cell_filtering = function(ref_rna, cell_type_col, filter_cell_type, min_cell_type_count){
  # 0. filter low count cell type
  low_count_cell_types = as.data.frame(table(ref_rna@meta.data[[cell_type_col]])) %>% 
                    setNames(c('cell_type','n')) %>% 
                    filter(n<min_cell_type_count) %>% pull(cell_type)

  # Filter
  print(filter_cell_type)
  print(length(low_count_cell_types)>0)
  if(filter_cell_type & length(low_count_cell_types)>0){
    print(str_glue("filtering out {toString(low_count_cell_types)}: cell count less than {min_cell_type_count}"))
    Idents(ref_rna) = cell_type_col
    ref_rna         = subset(ref_rna, idents = low_count_cell_types, invert = TRUE)
  }else if(!filter_cell_type & length(low_count_cell_types)>0){
    print(str_glue("Warning! Not filtering out {toString(low_count_cell_types)}: cell count less than {min_cell_type_count}"))
    print(str_glue("Might casue issue later"))
  }else{
    print(str_glue("Cell type counts all passed threshold {min_cell_type_count}"))
  }
  return(ref_rna)
}


# Make scRNA object ready for RCTD (Reference)
seurat_to_RCTD_ref = function(ref_rna, cell_type_col='cell_type', filter_cell_type =T, min_cell_type_count =25){
  # Check
  if(!cell_type_col %in% names(ref_rna[[]])){ 
    print(c("Availalbe names are: ", str_c(names(ref_rna[[]]),collapse = ', ') ))
    stop(str_glue('{cell_type_col} not found in the scRNA obj.'))
  }

  # Filtering
  ref_rna = RCTD_ref_cell_filtering(ref_rna=ref_rna, cell_type_col=cell_type_col, filter_cell_type =filter_cell_type, min_cell_type_count =min_cell_type_count)
  

  # 1. # meta_data  (with 3 columns, with headers “barcode”, “cluster”, and “nUMI”)
  meta_data = data.frame(barcode = colnames(ref_rna),
                         cluster = as.numeric(as.factor(ref_rna@meta.data[[cell_type_col]])) ,
                         nUMI    = ref_rna$nCount_RNA) %>% 
    filter(!is.na(cluster))
  
  # 2. cell_type_dict
  cell_type_dict = 
    data.frame(Cluster = as.numeric(as.factor(ref_rna@meta.data[[cell_type_col]])),
               Name    = ref_rna[[cell_type_col]]  ) %>% distinct %>% 
    setNames(c("Cluster","Name")) %>% 
    filter(!is.na(Cluster))
  
  
  # 3. dge: a Digital Gene Expression (DGE) (barcodes by gene counts) CSV file in the standard 10x format
  dge = GetAssayData(ref_rna) 
  
  # print(head(meta_data))
  # print(cell_type_dict)
  # print(dge[1:3,1:3])
  
  # Make reference seurat
  reference <- dgeToSeurat_direct(dge, meta_data, cell_type_dict)

  print('sc/snRNA Reference created successfully.')
  reference
}               


seurat_to_RCTD_ref_v2 = function(ref_rna, cell_type_col='cell_type', filter_cell_type =T, min_cell_type_count =25){
  # Updated 4/22/2021
  # 0. Check
  if(!cell_type_col %in% names(ref_rna[[]])){ 
    print(c("Availalbe names are: ", str_c(names(ref_rna[[]]),collapse = ', ') ))
    stop(str_glue('{cell_type_col} not found in the scRNA obj.'))
  }

  # 0. Filtering
  ref_rna = RCTD_ref_cell_filtering(ref_rna, cell_type_col=cell_type_col, filter_cell_type =filter_cell_type, min_cell_type_count =min_cell_type_count)

  # 1. Get data
  counts <- GetAssayData(object = ref_rna, slot = "counts") # load in counts matrix

  # 1.1 meta_data  (with 3 columns, with headers (barcodes, clusters, and nUMI)
  meta_data = data.frame(barcode = colnames(ref_rna),
                         cluster = as.factor(ref_rna@meta.data[['cell_type_Abbr']]) ,
                         nUMI    = colSums(counts)) %>% #ref_rna$nCount_RNA) %>% 
    filter(!is.na(cluster))

  # 1.2 Set named vectors
  cell_types = meta_data$cluster %>% setNames(meta_data$barcode)  # create cell_types named list
  nUMI       = meta_data$nUMI    %>% setNames(meta_data$barcode)  # create nUMI named list

  ### 2. Create the Reference object
  reference <- Reference(counts, cell_types, nUMI)

  print('sc/snRNA Reference created successfully.')
  reference
}




# Make ST object ready for RCTD
seurat_to_RCTD_st = function(st){
  # BeadLocationsForR.csv: a CSV file (with 3 columns, with headers “barcodes”, “xcoord”, and “ycoord”) containing the spatial locations of the pixels. 
  coords = lapply(Images(st), function(img) { 
    st@images[[img]]@coordinates[,c('imagecol','imagerow')]}) %>% 
    bind_rows() 
  
  # MappedDGEForR.csv: a DGE (gene counts by barcodes) CSV file. Represents raw counts at each pixel. 
  print('Getting DGE table')
  counts = GetAssayData(st_merge@assays$Spatial) # %>% as.data.frame %>% rownames_to_column('GENE')
  
  if(class(counts) != 'dgCMatrix'){
    message("Converting expression matirx to dgCMatrix")
    counts = as(as(counts, "matrix"), "dgCMatrix")
  }
  # Make st_merge obj
  print('Making ST object')
  puck = RCTD:::SpatialRNA(coords, counts)
  puck = RCTD:::restrict_puck(puck, colnames(puck@counts))
  print('ST quary for RCTD created successfully.')
  puck
}


get_RCTD_proportion = function(RCTD_obj, st_obj, add_to_st=T, assay_name='cell_type_RCTD'){
  if(missing(RCTD_obj)|missing(st_obj)){stop("Please provide both RCTD and ST object")}
  # result part
  results      = RCTD_obj@results
  # normalize the cell type proportions to sum to 1.
  norm_weights = sweep(results$weights, 1, rowSums(results$weights), '/')
  # Assign cell type percentage
  cell_weight_df = t(as.matrix(norm_weights))
  
  # Check if there's missing cell
  missing_cell  = setdiff(colnames(st_obj), rownames(norm_weights))
  if(length(missing_cell)>0){ # Missing cell
    # Create 0 df
    print(str_glue("Found {length(missing_cell)} missing spot(s), filling in with 0"))
    missing_cell_df = matrix(0, ncol = length(missing_cell), nrow = ncol(norm_weights)) %>% as.data.frame()
    rownames(missing_cell_df) = colnames(norm_weights)
    colnames(missing_cell_df) = missing_cell
    # Add to final cell weight df
    cell_weight_df = cbind(cell_weight_df, missing_cell_df)[,colnames(st_obj)]
  }else{
    cell_weight_df = cell_weight_df[,colnames(st_obj)]
  }
  
  # If add to ST or output as dataframe
  if(!add_to_st){ # output as dataframe
    print("Output as Dataframe mode.")
    return(cell_weight_df)
  }else{ # add to st_obj
    print(str_glue("Save to ST obj as {assay_name} assay"))
    st_obj[[assay_name]] = CreateAssayObject(cell_weight_df)  
    return(st_obj)
  }
}




# Internal function rewrite - createRCTD obj
# For debuging  - can use one from the package
# This version can handle case with "only 1 platform DEG" ---------=====
myrestrict_counts = function (puck, gene_list, UMI_thresh = 1, UMI_max = 20000)
{
  keep_loc = (puck@nUMI >= UMI_thresh) & (puck@nUMI <= UMI_max)
  new_matrix = puck@counts[gene_list, keep_loc] # filter counts
  if(length(gene_list)==1){ # only 1 gene left
    new_matrix = as(t(as.matrix(new_count)), 'dgCMatrix') # Manully convert back to sparce matrix format
    rownames(new_matrix) = gene_list
  }
  puck@counts = new_matrix
  if (length(puck@cell_labels) > 0)
    puck@cell_labels = puck@cell_labels[keep_loc]
  puck@nUMI = puck@nUMI[keep_loc]
  return(puck)
}

mycreate.RCTD =  function (spatialRNA, reference, max_cores = 8, test_mode = FALSE,
                           gene_cutoff = 0.000125, fc_cutoff = 0.5, gene_cutoff_reg = 2e-04,
                           fc_cutoff_reg = 0.75, UMI_min = 100, UMI_max = 2e+05, class_df = NULL,
                           CELL_MIN_INSTANCE = 25)
{
  config <- list(gene_cutoff = gene_cutoff, fc_cutoff = fc_cutoff,
                 gene_cutoff_reg = gene_cutoff_reg, fc_cutoff_reg = fc_cutoff_reg,
                 UMI_min = UMI_min, max_cores = max_cores, N_epoch = 8,
                 N_X = 50000, K_val = 100, N_fit = 1000, N_epoch_bulk = 30,
                 MIN_CHANGE_BULK = 1e-04, UMI_max = UMI_max, MIN_OBS = 3)
  if (test_mode)
    config <- list(gene_cutoff = 0.00125, fc_cutoff = 0.5,
                   gene_curoff_reg = 0.002, fc_cutoff_reg = 0.75, UMI_min = 1000,
                   N_epoch = 1, N_X = 50000, K_val = 100, N_fit = 50,
                   N_epoch_bulk = 4, MIN_CHANGE_BULK = 1, UMI_max = 2e+05,
                   MIN_OBS = 3, max_cores = 1)
  cell_type_info <- list(info = RCTD:::process_cell_type_info(reference,
                                                              CELL_MIN = CELL_MIN_INSTANCE), renorm = NULL)
  puck = myrestrict_counts(spatialRNA, rownames(spatialRNA@counts),
                           UMI_thresh = config$UMI_min, UMI_max = config$UMI_max)
  print('Here0')
  
  # DEGs
  print("create.RCTD: getting regression differentially expressed genes: ")
  gene_list_reg = RCTD:::get_de_genes(cell_type_info$info, puck, 
                                      fc_thresh = config$fc_cutoff_reg,
                                      expr_thresh = config$gene_cutoff_reg, MIN_OBS = config$MIN_OBS)
  print(gene_list_reg %>% head)
  if (length(gene_list_reg) == 0)
    stop("create.RCTD: Error: 0 regression differentially expressed genes found")
  
  # Platform effect
  print("create.RCTD: getting platform effect normalization differentially expressed genes: ")
  print('HERE1')
  gene_list_bulk = RCTD:::get_de_genes(cell_type_info$info, puck,
                                       fc_thresh = config$fc_cutoff, expr_thresh = config$gene_cutoff,
                                       MIN_OBS = config$MIN_OBS)
  # gene_list_bulk = c('MIR1302-2HG',gene_list_bulk) # Adding a muck one
  print(gene_list_bulk %>% head)
  
  if (length(gene_list_bulk) == 0)
    stop("create.RCTD: Error: 0 bulk differentially expressed genes found")
  
  puck = myrestrict_counts(puck, gene_list_bulk, UMI_thresh = config$UMI_min,
                           UMI_max = config$UMI_max)
  
  puck = RCTD:::restrict_puck(puck, colnames(puck@counts))
  if (is.null(class_df))
    class_df <- data.frame(cell_type_info$info[[2]], row.names = cell_type_info$info[[2]])
  colnames(class_df)[1] = "class"
  internal_vars <- list(gene_list_reg = gene_list_reg, gene_list_bulk = gene_list_bulk,
                        proportions = NULL, class_df = class_df)
  
  print("RCTD created successfully.")
  new("RCTD", spatialRNA = puck, reference = reference, config = config,
      cell_type_info = cell_type_info, internal_vars = internal_vars)
}
# ------------==========