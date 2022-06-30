# This script make 1 page summary of the output
library(tidyverse)
library(Seurat)
library(patchwork)
library(grid)


make_ST_summary_plot = function(folder_path, out_path, plot_title, replace_exist = F, seurat_out_path){
  ptm = Sys.time()
  # process path
  if(!'outs' %in% list.files(folder_path) ) stop('outs folder not found. Please double check path')
  if(str_sub(folder_path,-1) == "/") folder_path = str_sub(folder_path, 1, -2) # remove ending "/"
  if(missing(out_path)){out_path = folder_path
  }else{
    if(str_sub(out_path,-1) == "/") out_path = str_sub(folder_path, 1, -2) # remove ending "/"
  }

  # Load metrics first to check files
  metrics  = data.table::fread(str_glue('{folder_path}/outs/metrics_summary.csv'))
  # set value
  if(missing(plot_title)) plot_title = metrics$`Sample ID`
  output_file = str_glue('SummaryPlot_{metrics$`Sample ID`}.pdf')
  # Check if output exist
  if(!replace_exist){
    if(output_file %in% list.files(out_path)){
      message(str_glue('Figure {output_file} found in {out_path}. Skip.'))
      return()
    }
  }
  
  # Loading
  message('Loading files')
  img_g    = rasterGrob(png::readPNG(str_glue("{folder_path}/outs/spatial/tissue_hires_image.png")), interpolate=TRUE)
  obj      = Seurat::Load10X_Spatial(str_glue('{folder_path}/outs'))
  
  
  # Seurat processing
  message('Procesing Seurat')
  obj <- SCTransform(obj, assay = "Spatial", verbose = FALSE) %>%
    RunPCA( assay = "SCT", verbose = FALSE) %>%
    FindNeighbors( reduction = "pca", dims = 1:30) %>%
    FindClusters( verbose = FALSE) %>%
    RunUMAP( reduction = "pca", dims = 1:30)
  
  p1 <- VlnPlot(obj, features = "nCount_Spatial", pt.size = 0.1) + labs(x='',y='nCount_Spatial',title = '') + theme(axis.title.y = element_text(size=15, face = 'bold', vjust = 0.5))
  p2 <- SpatialFeaturePlot(obj, features = "nCount_Spatial", pt.size.factor = 1.8,stroke=NA)
  p3 <- DimPlot(obj, reduction = "umap", label = TRUE) + labs(x='',y='UMAP',title = '') + theme(axis.title.y = element_text(size=15, face = 'bold', vjust = 0.5),
                                                                                                axis.text = element_blank(),
                                                                                                axis.ticks = element_blank())
  p4 <- SpatialDimPlot(obj, stroke=NA, label=T, label.box = F, label.size = 4, pt.size.factor = 1.8)
  design1 = "CD
             AB"
  p_all = wrap_plots(C=p3, D=p4, A=p1, B=p2, design = design1) & NoLegend()
  
  
  subtitle = str_glue("Number of Spots Under Tissue: {metrics$`Number of Spots Under Tissue`} | Mean Reads per Spot: {round(metrics$`Mean Reads per Spot`,0)} | Median Genes per Spot: {metrics$`Median Genes per Spot`}")
  
  # Plot result
  
  pdf(str_glue('{out_path}/{output_file}'), width = 13, height = 8)
  ## V1
  p = (p_all | wrap_elements(img_g) ) + plot_annotation(
    title = plot_title,
    subtitle = subtitle,
    theme = theme(plot.title = element_text(size = 30, face = 'bold'),
                  plot.subtitle = element_text(size = 15))
  )
  
  #@@ V2
  # design = "AB
  #           CC"
  # p_he = wrap_elements(img_g)
  # p_side = wrap_plots(A = p1, B=p3, C=p_he, design = design) & NoLegend()
  # p =((p2 | p4 | p_side) & NoLegend() )+ plot_annotation(
  #   title = plot_title,
  #   subtitle = subtitle,
  #   theme = theme(plot.title = element_text(size = 30, face = 'bold'),
  #                 plot.subtitle = element_text(size = 15))
  # )
  #@@ V3
  # design = 'AAAABB'
  # p = wrap_plots(A = p_all, B= wrap_elements(img_g), design = design) + plot_annotation(
  #   title = plot_title,
  #   subtitle = subtitle,
  #   theme = theme(plot.title = element_text(size = 30, face = 'bold'),
  #                 plot.subtitle = element_text(size = 15))
  # )
  print(p)
  dev.off()
  if(!missing(seurat_out_path)){
    message(str_glue("Saving processing Seurat obj. Outputed to {seurat_out_path}/{metrics$`Sample ID`}.rds"))
    saveRDS(obj, str_glue('{seurat_out_path}/{metrics$`Sample ID`}.rds'))
    
  }
  Sys.time()-ptm;ptm = Sys.time()
  message(str_glue("Completed. Outputed to {out_path}/SummaryPlot_{metrics$`Sample ID`}.pdf"))
}

# EXAMPLE: 
# make_ST_summary_plot('~/Box/Ding_Lab/Projects_Current/Spatial_transcriptomics/processed_data/NMK113F/NMK113F/')
