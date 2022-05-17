# This is modified from Ruiyang's dot_plot function
# V2 
library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

# Functions
# Helper
helper_cluster_expression_plot = function(to_plot, exprs.to.plot, dot.color.range, facet.order, y.axis.order, rotate_title){
	## plotting
  if (!is.null(facet.order)){to_plot$cell_type <- to_plot$cell_type %>% factor(levels=facet.order)}
  if (!is.null(y.axis.order)){to_plot$id <- to_plot$id %>% factor(levels=y.axis.order)}
  p <- ggplot(to_plot)+
    geom_point(aes_string(x="features.plot",y="id",size="pct.exp",color=exprs.to.plot),shape=16)+scale_color_gradientn(colors=dot.color.range)+
    facet_grid(.~cell_type,scale="free",space="free")+
    theme(panel.background=element_blank(),panel.border=element_rect(fill=NA,color="black"))+
    theme(axis.text.x=element_text(angle=90))+
    theme(strip.background=element_rect(fill=NA),strip.text.x=element_text(size=10,face="bold"))+
    theme(axis.text.x=element_text(vjust=0.5,hjust=1))
  if(rotate_title){p = p + theme(strip.text=element_text(angle=90, vjust=0.5, hjust=0))}
  return(p)
}

helper_check_celltype_marker = function(celltype_markers){
  ## celltype_markers should be a data frame where column names include "gene_symbol" and "cell_type"
  if(is.null(celltype_markers)){
    stop('Need cell type marker dataframe include "gene_symbol" and "cell_type" column')
  }
  if(!'gene_symbol'%in%colnames(celltype_markers)|| !'cell_type'%in%colnames(celltype_markers)){
    stop('Celltype_markers should be a data frame where column names include "gene_symbol" and "cell_type"')
  }
}

helper_filter_cell_marker = function(obj, celltype_markers, min_marker=1){
  # Filter cell type with enough markers (default =1, set 0 to ignore)
  celltype_markers = filter(celltype_markers, gene_symbol%in%rownames(obj))
  message(str_glue('Removing cell type with less than {min_marker} marker(s)'))
  celltype_markers = filter(celltype_markers, 
                            cell_type%in%filter(count(celltype_markers, cell_type), 
                                                n>min_marker)$cell_type)
  celltype_markers
}

# Main functions -----------------
cluster_DotPlot <- function(Seurat_obj,celltype_markers=NULL,group.by="seurat_clusters",exprs.to.plot="avg.exp.scaled",
											filter_data=F, pct.filter=25,avg.exp.filter=1, min_marker=1,
                                            dot.color.range=brewer.pal(9,"YlOrRd")[3:9],facet.order=NULL,y.axis.order=NULL,
                                            rotate_title=F, add_summary=T){
  helper_check_celltype_marker(celltype_markers)
  ## exprs.to.plot: whether to choose avg.exp or avg.exp.scaled for plotting (shown by color)
  celltype_markers = helper_filter_cell_marker(Seurat_obj, celltype_markers, min_marker)
  
  # Make plot df
  p <- DotPlot(Seurat_obj,features=celltype_markers$gene_symbol %>% as.character %>% unique,group.by=group.by)
  to_plot <- p$data
  to_plot <- merge(to_plot,celltype_markers,by.x="features.plot",by.y="gene_symbol",all.x=TRUE)

  
  ## annotation if some genes are markers for multiple cell types
  duplicated_index <- grep(TRUE,duplicated(celltype_markers$features.plot %>% as.character))
  if (length(duplicated_index)>0){
    to_plot$features.plot[duplicated_index] <- paste0(to_plot$features.plot,"_",to_plot$cell_type)[duplicated_index]
  }
  
  # Filter 
  if(filter_data){
  	message(str_glue('Filtering data by pct.filter={pct.filter}, avg.exp.filter = {avg.exp.filter}'))
    genes_to_keep  = to_plot %>% group_by(features.plot) %>%
      summarize(max.pct=max(pct.exp),max.exp=max(avg.exp)) %>%
      filter(max.pct>=pct.filter & max.exp>=avg.exp.filter) %>% .$features.plot
  	message(str_glue("Removing {length(setdiff(to_plot$features.plot, genes_to_keep))} genes. Plotting {length(genes_to_keep)} genes."))
  	to_plot = to_plot %>% filter(features.plot %in% genes_to_keep)
  }

  
  # Make summary dot plot
  summary_dot_df = to_plot %>% group_by(cell_type, id) %>% summarize(avg.exp = mean(avg.exp), pct.exp= mean(pct.exp), avg.exp.scaled= mean(avg.exp.scaled))

  ## plotting
  p = helper_cluster_expression_plot(to_plot, exprs.to.plot, dot.color.range, facet.order, y.axis.order, rotate_title)
  
  # Add summary dotplot
  if(add_summary){
    p_summary = ggplot(summary_dot_df, aes(x = cell_type, y=id, size = pct.exp, color = avg.exp.scaled,)) + 
      geom_point() + 
      scale_color_gradientn(colors=dot.color.range)+
      theme_bw() + 
      theme(axis.text.x=element_text(angle=-45, vjust = 0.5, hjust = 0, face='bold')) + labs (x="")
    p = p+p_summary+plot_layout(widths = c(9, 1))
  }
  return(p)
}






# Cell score function
getCellTypeMarkerAssay = function(obj, celltype_marker=NULL, new_assay_name="CellTypeScore", min_marker = 1){
  DefaultAssay(obj)='SCT'
  # cell type marker process
  helper_check_celltype_marker(celltype_marker)
  celltype_marker = filter(celltype_marker, gene_symbol %in% rownames(obj))
  celltype_marker = helper_filter_cell_marker(obj, celltype_marker, min_marker)
  
  # Normalize expression of markers
  expression = GetAssayData(obj)[celltype_marker$gene_symbol, ] %>% as.matrix
  expression = t(scale(t(expression)))
  
  # Calculate average normalize expression per cell_type
  tmp=lapply(unique(celltype_marker$cell_type), function(type){
    tmp_marker = filter(celltype_marker, cell_type==type)$gene_symbol
    tmp_exp = filter(expression %>% as.data.frame %>% rownames_to_column('gene'), gene %in% tmp_marker) 
    tmp_exp = tmp_exp %>% select(-gene)
    if(ncol(tmp_exp)>1){
      colMeans(as.matrix(tmp_exp))
    }else{ # Only 1 gene in this cell type. Retun as is
      tmp_exp 
    }
  })
  
  # Generate cell type score dataframe
  score_df = bind_rows(tmp) %>% as.data.frame()
  rownames(score_df) = unique(celltype_marker$cell_type)
  
  # save to one of the Seurat assay
  obj[[new_assay_name]] = CreateAssayObject(counts = score_df)
  print(str_glue('Output saved in Seurat obj as a new assay named - {new_assay_name}'))
  return(obj)
}

