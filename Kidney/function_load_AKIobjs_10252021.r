# Functions to load AKI object
# updated 10/25/2021
library(tidyverse)
library(qs)
library(Seurat)

load_nmk_st = function(){
        ## ST files
        project_folder = "/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/NMK/NMK_integration_v5"
        st_path    = 'obj/NMK_merge_v5_11042021.qs'

        message(str_glue('This version load ST from the following:
                        obj : {st_path}'))
        # Load
        obj = qread(st_path)
        
        return(obj)
    }


loadakiobj = function(type = c('rna','atac','st')){
    type = match.arg(type)
    switch(type,
           rna = loadakisn(),
           st = loadakiST(),
           atac = loadakiatac(),
          )
}

loadakisn = function(){
        print('Loading snRNA: Current version updated on 10/25/2021')
        snaki  = readRDS("/diskmnt/Projects/Users/rliu/Projects/normal_mouse_kidney_development/snRNA/Seurat_merging/AKI/AKI_W12X2_3d_8d/analysis/doublet_removal/results/integration_AKI_W12X2_3d_8d_doublets_removed.rds")
        snmeta = read.table('/diskmnt/Projects/Users/rliu/Projects/normal_mouse_kidney_development/snRNA/Seurat_merging/AKI/AKI_W12X2_3d_8d/analysis/analysis_after_doublet_removal/meta.data/results/meta.data_20210810_v2.txt', header = T)
        snaki@meta.data = snmeta

        ## Add meta data
        snaki@meta.data = snaki@meta.data %>% 
                        mutate(Sex_cell_group = str_c(Sex, cell_type_Abbr, treatment_group, sep = '__'))
        return(snaki)
    }

loadakiST = function(){
        print('Loading ST: Current version updated on 09/20/2021')
        ## ST files
        aki_project_folder = "/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/NMK/NMK_integration_v4"
        st_path    = str_glue('{aki_project_folder}/obj/NMK_AKI_merge_v4_04142021.qs')
        meta_path  = 'table/Meta_AKI_v1_RCTD_doubletRM_matched_09192021.rds'
        assay_path = 'table/Assay_AKI_v1_RCTD_doubletRM_matched_09192021.rds'

        message(str_glue('This version load ST from the following:
                        AKI : {st_path}
                        meta : {meta_path}
                        assay : {assay_path}'))
        # Load
        AKI = qread(st_path)
        ## Assay
        assay_obj = readRDS(assay_path)
        AKI$RCTD_doubletRM_v1 = assay_obj
        ## meta
        meta_new = readRDS(meta_path)
        AKI@meta.data = meta_new[rownames(AKI@meta.data), ]
    
        return(AKI)
    }

loadakiatac = function(){
    library(Signac)
    atac = readRDS('/diskmnt/Projects/Users/rliu/Projects/normal_mouse_kidney_development/snATAC/Merging/AKI/AKI_W12X2_3d_8d/with_new_clusters_pc_2_30/AKI_W12X2_3d_8d_snATAC_Merged.rds')
    # Add info
    # ATAC
    atac@meta.data = atac@meta.data %>% mutate(Sex = case_when(str_detect(dataset, 'F-|F1-')~'F',
                                              str_detect(dataset, 'M-|M1-')~'M'),
                              Day = case_when(str_detect(dataset, 'AKI8')~'8day',
                                              str_detect(dataset, 'AKI3')~'3day',
                                              T~'Control'),
                              Day = factor(Day, levels = c('Control', '3day', '8day')),
                            Sex_day = str_c(Sex, Day, sep='_'),
                            Sex_day_cell = str_c(Sex, Day, cell_type_Abbr_scRNA, sep='__'),
                                    )
    return(atac)
    }

