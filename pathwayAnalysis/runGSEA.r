message("conda activate clusterProfiler")
# Convert Symbol to ENTREZID 
library(clusterProfiler)
library(msigdbr)
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)
message("Using organism: ", organism)
message("Available functions: gene2entrezid, genevectorName2entrezid, genevectorName2entrezid")
gene2entrezid = function(gene, organ_database = "org.Hs.eg.db", na.rm =F){
    eg_df = bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb=organ_database)
    converted_id = eg_df$ENTREZID %>% setNames(eg_df$SYMBOL) 
    if(na.rm) converted_id = converted_id %>% .[!is.na(.)]
    return(converted_id)
}
genevectorName2entrezid = function(fc_vec_genelist,  organ_database = "org.Hs.eg.db", na.rm =T){
    gene_symbols = names(fc_vec_genelist)
    entrezid_vec = gene2entrezid(gene_symbols)
    gene_list_w_entrezid = fc_vec_genelist %>% setNames(entrezid_vec[gene_symbols])
    if(na.rm) gene_list_w_entrezid = gene_list_w_entrezid %>% .[!is.na(names(.))]
    return(gene_list_w_entrezid)
}

# Wrapper From DEG table to GSEA result
message("Wrapper functions: FindMarkerTable2Genelist, FindAllMarkerTable2Genelist")
# Convert from FindMarker table to gene vector.
# assume rownames already convert to column called "gene"
FindMarkerTable2Genelist = function(deg_df, gene_colname = 'gene', log2FC_colname = 'avg_log2FC'){
    # Prepare input format and remove cluster with too few DEGs
    # avg_log2FC as values, gene names as vector
    # Sort value decreasing
    deg_df %>% {setNames(.[[log2FC_colname]], .[[gene_colname]])} %>% sort(., decreasing = T)
}

# Conver from FindAllMarker table to list of gene vector 
# Assume the "cluster" column have each identity
FindAllMarkerTable2Genelist = function(deg_df, split_column = 'cluster', gene_colname = 'gene', log2FC_colname = 'avg_log2FC', n_min_deg = 5){
    # Prepare input format and remove cluster with too few DEGs
    message("Minimum number of DEG:", n_min_deg)
    deg_list= deg_df %>% 
        split(.[[split_column]]) %>% 
        # avg_log2FC as values, gene names as vector, Sort value decreasing
        map(~setNames(.x[[log2FC_colname]], .x[[gene_colname]]) %>% sort(decreasing = T)) 

    # Remove group with less than n deg
    deg_list %>% keep(~ length(.) > n_min_deg)
}   

## ------------------ GSEA ------------------
# Wrapper for Hallmark
message("Wrapper functions: FindMarkerTable2GSEAresult, FindAllMarkerTable2GSEAresult")
message("Note, default, the wrapper function uses geneset: geneset_df_HALMARK. If need replacement, please use the argument geneset_use")
message("See this page section 12.3 of how to get the dataset: http://yulab-smu.top/biomedical-knowledge-mining-book/universal-api.html")
geneset_df_HALMARK <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  #dplyr::select(gs_name, entrez_gene) %>% 
  dplyr::select(gs_name, human_gene_symbol) %>% # New, use human gene symbol
    # Remove HALLMARK
    mutate(gs_name = str_remove(gs_name, "HALLMARK_"))


ShowMsigDBHumanCategories = function(){
    message("H: hallmark gene sets
    C1: positional gene sets
    C2: curated gene sets
    C3: motif gene sets
    C4: computational gene sets
    C5: GO gene sets
    C6: oncogenic signatures
    C7: immunologic signatures")
}

LoadMSigDBHuman = function(){
    message("Loading All Human MsigDB geneset into geneset_msig_ALLHUMAN_df, geneset_msig_ALLHUMAN_use_df, geneset_msig_ALLHUMAN_list")
    geneset_msig_ALLHUMAN_df <<- msigdbr(species = "Homo sapiens") 
    geneset_msig_ALLHUMAN_list <<- geneset_msig_ALLHUMAN_df %>% split(f = .$gs_cat) %>% map(., ~dplyr::select(., gs_name, human_gene_symbol))
    ShowMsigDBHumanCategories()
    geneset_msig_ALLHUMAN_use_df <<- geneset_msig_ALLHUMAN_df %>% dplyr::select(gs_name, human_gene_symbol) 
    message("Use geneset_use = geneset_msig_ALLHUMAN_use_df in FindMarkerTable2GSEAresult")
}

FindMarkerTable2GSEAresult = function(deg_df, gene_colname = 'gene', log2FC_colname = 'avg_log2FC', geneset_use = geneset_df_HALMARK, scoreType = "std"){
    # Preflight checks ..
    gene_use = FindMarkerTable2Genelist(deg_df)
    message("Length of gene:", length(unique(deg_df[[gene_colname]])))
    #gene_vec_converted = genevectorName2entrezid(gene_use)
    # New, uses Huamn Gene symbol as reference now
    gene_vec_converted = gene_use 
    message("Length of converted:", length(gene_vec_converted))

    # first check if any gene in the geneset ID
    #geneset_overlap = intersect(geneset_use$entrez_gene, names(gene_vec_converted))
    geneset_overlap = intersect(geneset_use$human_gene_symbol, names(gene_vec_converted))
    message("Length of intersect:", length(geneset_overlap))
    if(length(geneset_overlap) == 0){message("No intersect found. Skip"); return(NULL)}
    
    # Run GSEA
    GSEA(gene_vec_converted, TERM2GENE = geneset_use, scoreType = scoreType)
}
# Run genesets from MsigDB
FindMarkerTable2GSEA_MsigDB = function(deg_df, genesets_include = c('H','C1','C2','C3','C4','C5','C6'), ...){
    ShowMsigDBHumanCategories()
    message("Using following MSigDB genesets:", toString(genesets_include))

    # Load database if not yet
    if(!exists("geneset_msig_ALLHUMAN_list" )) LoadMSigDBHuman()

    geneset_msig_list_use = geneset_msig_ALLHUMAN_list[genesets_include]
    imap(geneset_msig_list_use, function(geneset_df, cat_name){
        message("Running: ", cat_name)
        FindMarkerTable2GSEAresult(deg_df, geneset_use = geneset_df, ...)
    })
}




FindAllMarkerTable2GSEAresult = function(deg_df, split_column = 'cluster', 
    gene_colname = 'gene', log2FC_colname = 'avg_log2FC', 
    n_min_deg = 5, remove_empty_results = T,
    geneset_use = geneset_df_HALMARK,
    scoreType = "std"
    ){
    find_all_marker_gene_list = FindAllMarkerTable2Genelist(deg_df, split_column, gene_colname, log2FC_colname, n_min_deg)
    result_list = imap(find_all_marker_gene_list, function(gene_use, group_name){
        message("Running GSEA for group:", group_name)
        message("Length of gene:", length(gene_use))
        #gene_vec_converted = genevectorName2entrezid(gene_use)
        # New, uses Huamn Gene symbol as reference now
        gene_vec_converted = gene_use 

        message("Length of converted:", length(gene_vec_converted))
        # first check if any gene in the geneset IDr
        #geneset_overlap = intersect(geneset_use$entrez_gene, names(gene_vec_converted))
        geneset_overlap = intersect(geneset_use$human_gene_symbol, names(gene_vec_converted))
        message("Length of intersect:", length(geneset_overlap))
        if(length(geneset_overlap) == 0){message("No intersect found. Skip"); return(NULL)}
        print(head(gene_vec_converted))
        GSEA(gene_vec_converted, TERM2GENE = geneset_use, scoreType = scoreType)
    })

    # Clean up result
    if(remove_empty_results){
        message("Remove empty results")
        result_list = imap(result_list, function(gsea_result, clusterID){
            if(is.null(gsea_result)) return(NULL)
            if(nrow(gsea_result@result) == 0){
                # No GSEA result
                message("No GSEA result in ", clusterID)
                return(NULL)
            }
            message("Have GSEA result in ", clusterID)
            return(gsea_result)
        }) %>% discard(is.null)
    }  
    return(result_list)
}


# Run genesets from MsigDB
FindAllMarkerTable2GSEAresult_MsigDB = function(deg_df, genesets_include = c('H','C1','C2','C3','C4','C5','C6'), ...){
    ShowMsigDBHumanCategories()
    message("Using following MSigDB genesets:", toString(genesets_include))

    # Load database if not yet
    if(!exists("geneset_msig_ALLHUMAN_list" )) LoadMSigDBHuman()

    geneset_msig_list_use = geneset_msig_ALLHUMAN_list[genesets_include]
    imap(geneset_msig_list_use, function(geneset_df, cat_name){
        message("Running: ", cat_name)
        FindAllMarkerTable2GSEAresult(deg_df, geneset_use = geneset_df, ...)
    })
}


## ---- Utility functions ----
# Make dotplot compare Tumor and TME
combineGSEAresultslist = function(gsea_result_list){
    gsea_result_list %>% imap(~as.data.frame(.x) %>% rownames_to_column('geneset') %>% mutate(group = .y) 
    )  %>% bind_rows() 
}

ReorderByHCluster = function(data_use, ident_column, groupby_column, values_column){
    # 1. convert to wide
    value_mtx = data_use %>% 
        dplyr::select(all_of(c(ident_column,groupby_column, values_column))) %>%
        pivot_wider(id_cols = all_of(ident_column), names_from = all_of(groupby_column), values_from = all_of(values_column), values_fill=0) %>% 
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

## ---- Plot functions ----
MakeGeneSetDotplot = function(gsearesult_list, plt_title = "GSEA result", subtitle = "Geneset", reorder_geneset = T, reorder_group = T, ...){
    combined_df = combineGSEAresultslist(gsearesult_list)
    # Make dotplot x = group, y = geneset, color = enrichmentScore, size = -log10(qvalue)
    # Reorder geneset/group by hclust
    if(reorder_geneset) combined_df = ReorderByHCluster(data_use = combined_df, 'group', 'geneset', 'NES') 
    if(reorder_group) combined_df = ReorderByHCluster(data_use = combined_df, "geneset", 'group', 'NES')
        
    ggplot(combined_df, aes(x = group, y = geneset, color = NES, size = -log10(qvalue))) + geom_point() + 
        scale_color_gradient2(low = "#2222EE", mid = 'gray80', high = "#FF4444") + 
        theme_bw() + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
        #facet_wrap(~group, scales = 'free_y') + 
        labs(title = plt_title, subtitle = subtitle, color = "Normalized\nEnrichment\nScore", size = '-log10(qvalue)') 
        #theme(legend.position = 'none')
}