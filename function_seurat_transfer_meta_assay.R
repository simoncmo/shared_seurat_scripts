# Add Meta and ASSAY after transfer
# For inidividual object and Merge one


#### EXAMPLE USE ########
# A. individual obj list
# obj_list_cell = AddMetaAssayToObjList(obj_list, predict_cellTypeAbbr_list, assay_name = 'cell_type_Abbr') %>% 
#      AddMetaAssayToObjList(obj_list = ., predict_cellTypeAfterTransfer_list, assay_name = 'cell_type_after_transfer') 

# B. Merged list
# merged_new = AddMetaAssayListToMergedObj(merged_obj, predict_cellTypeAbbr_list, assay_name ='cell_type_Abbr') %>% 
#     AddMetaAssayListToMergedObj(predict_cellTypeAfterTransfer_list, assay_name ='cell_type_after_transfer')
########################
### Preiction to OBJ FUNCTIONs
# - Clean Cell name
# - Get max cell
# - Create obj

# - MERGE ONLY
#    - Add sample name
#    - Merge matrix

## WRAPPER ################
AddMetaAssayToObjList = function(obj_list, predition_list, assay_name='new'){
    map2(obj_list, predition_list, function(obj, predict){
        new_assay = predict %>% CleanCellName %>% CreateCleanAssay
        new_meta  = predict %>% CleanCellName %>% GetMaxCell(new_name = assay_name)

        obj = AddMultipleMeta(obj, new_meta)
        obj[[str_glue('{assay_name}_assay')]] = new_assay
        obj
    })
}

AddMetaAssayListToMergedObj = function(merged_obj, predition_list, assay_name='new'){
    # Merge list
    merged_df = predition_list %>% CleanCellName %>% MakeMergedPredictionMatrix 
    # Meta
    merged_obj = AddMultipleMeta(merged_obj, merged_df %>% GetMaxCell(new_name = assay_name))
    # Asay
    merged_obj[[str_glue('{assay_name}_assay')]] = CreateCleanAssay(clean_prediction = merged_df)
    return(merged_obj)
}



##########################
# ------ Individual STEPS
##########################
CleanCellName = function(prediction){
    if(inherits(prediction, 'list')) map(prediction, ~.x %>% setNames(str_remove(names(.), 'prediction.score.')) )
    else prediction %>% setNames(str_remove(names(.), 'prediction.score.')) 
}

GetMaxCell = function(clean_prediction, new_name){
    result = clean_prediction[,c('predicted.id','max')] 
    if(!missing(new_name)) result = result %>% setNames(c(new_name, str_glue('{new_name}_score')))
    result
}

CreateCleanAssay = function(clean_prediction){
    clean_prediction[,colnames(clean_prediction)!='predicted.id'] %>% 
        setNames(str_replace(names(.), '_','-')) %>% 
        t %>% 
        CreateAssayObject
}

AddMultipleMeta = function(obj, ...){
    # Meta data process
    meta_merged = bind_cols(list(...)) # meta to add to obj
    row_order   = obj@meta.data %>% rownames
    if(length(setdiff(row_order, rownames(meta_merged) ))!=0) stop('rowname mismatch with obj. Please double check')

    AddMetaData(object = obj, metadata = meta_merged[row_order, ])
}

## FOR MERGE
MakeMergedPredictionMatrix = function(prediction_list, sample_names){
    # Fix name and combine
    if(!missing(sample_names)) names(preidction_list) = sample_names
    imap(prediction_list, function(predict_mtx, sample){
        rownames(predict_mtx) = str_c(sample, rownames(predict_mtx), sep='_')
        predict_mtx
    }) %>% bind_rows()
}