# Plot RCTD output
library(Seurat)
library(qs)
library(tidyverse)
library(patchwork)
suppressPackageStartupMessages(library("argparse"))

project_folder = "/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/NMK/NMKST_v6"
setwd(project_folder)
source('NMKST_v6/script/Utility/Kidney_colorset_v4_12132021.R')
source('NMKST_v6/script/PLOT/function_Kidney_overview_plot_122021.r')
source('/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/Shared_resource/script/function_seurat_pipeline_v1.R')
source('/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/Shared_resource/script/function_seurat_helper.R')
source('/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/Shared_resource/script/function_spatial_plot_helper.R')
source('/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/Shared_resource/script/Kidney/function_Kidney_STDim_v5.R')

# Define col and assay names
process_date   = format(Sys.Date(), "%m%d%Y")

# Color
source("NMKST_v6/script/Utility/Kidney_colorset_v4_12132021.R")
source("NMKST_v6/script/Utility/color_fun.R")

###### Load Arg ####
###############################################################
# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("-s", "--mergedst", type="character", default=NULL, 
                    help="Merged ST object from RCTD")
parser$add_argument("-b", "--basepath", type="character", default=NULL, 
                    help="Base Path of Merged ST object from RCTD")
parser$add_argument("-c", "--cellcol", type="character", default=NULL,
                    help="Cell col name for the Plot")
parser$add_argument("-w", "--pdfwidth", type="integer", default=20, 
                    help="PDF output width")
parser$add_argument("-t", "--pdfheight", type="integer", default=7, 
                    help="PDF output height")
parser$add_argument("-o", "--out", type="character", default="STDistributionPlot", 
                    help="Figure output path [default %(default)s]")

# Get arguments
opt <- parser$parse_args()

## Process ARG
cell_col_name = str_c('RCTD',opt$cellcol, 'name', sep='_')
path_out = str_glue('{opt$basepath}/{opt$out}')


## LOAD
###############################################################
stnormal = qread(str_glue("{opt$basepath}/{opt$mergedst}"))
meta_stnormal = stnormal@meta.data

## Set up 
###############################################################
## Extend color
col_cell_all = ExtendColorVector(col_cell_type, stnormal@meta.data[[cell_col_name]])

## Order
order_age = c('E165','P0','W1','W2','W3','W12','W52','W92','W113')
stnormal@meta.data = stnormal@meta.data %>% mutate(Age = factor(Age, levels = order_age))


## PLOT
###############################################################
message('Making plot ..')
p = make_normal_kidney_plots(obj = stnormal, 
                             variables = cell_col_name, 
                             col_cell_all = col_cell_all, 
                             plt_mode ='ident')

## Save
# Outputs:
message(str_glue("Saving to {path_out}/STDistribution_{cell_col_name}.pdf"))
dir.create(str_glue('{path_out}'), recursive = T)
pdf(str_glue("{path_out}/STDistribution_{cell_col_name}.pdf"), width = opt$pdfwidth, height = opt$pdfheight)
p
dev.off()
