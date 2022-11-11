# SETUP
suppressPackageStartupMessages(library(tidyverse))
#suppressPackageStartupMessages(library(qs))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(spacexr))
suppressPackageStartupMessages(library(optparse))

source('/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/Shared_resource/script_git/Spacexr/function_spacexr_helper_v2.R')

# tiktoc
today_date = Sys.Date() %>% format('%m%d%Y')
ptm = Sys.time()

# Parse input
option_list = list(
    make_option(c("--sn"), type="character", default=NULL, 
              help="snRNA obj path", metavar="character"),
    make_option(c("--snmeta"), type="character", default=NULL, 
              help="seurat obj  meta path", metavar="character"),
    make_option(c("--st"), type="character", default=NULL, 
              help="ST obj path", metavar="character"),
    make_option(c("--sntransfer"), type="character", default=NULL, 
              help="snRNA column name to transfer", metavar="character"),
    make_option(c("--transfername"), type="character", default='transfer', 
              help="Name for this transfer", metavar="character"),
    make_option(c("-o", "--out"), type="character", default="RCTDresult", 
              help="output path [default= %default]", metavar="character"),
    make_option(c("--sample"), type="character", default=NULL, 
              help="Sample name", metavar="character"),
    make_option(c("--cellthreshold"), type="integer", default=25,
	      help="Cell count threshold for filtering")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# function obj loader
obj_loader = function(obj_path){
    message(str_glue('loading object from {obj_path}..'))
    READ_FUN = switch(str_extract(obj_path,'[[.]]rds|[[.]]RDS|[[.]]qs'), .rds=readRDS, .RDS=readRDS, .qs=qread)
    if(is.null(READ_FUN)){stop('No supported file type detected')}
    return(READ_FUN(obj_path))
}

# Load
sn = obj_loader(opt$sn)
ST = obj_loader(opt$st)

# Meta for snRNA
if(!is.null(opt$snmeta)){
    message('Adding meta data to snRNA obj')
    snmeta = read.table(opt$snmeta, header=T)
	snmeta = snmeta[rownames(sn@meta.data), ]
	sn@meta.data = snmeta[rownames(sn@meta.data), ]
    message(str_glue('Column after adding meta = {toString(names(sn@meta.data))}'))
}else{
	snmeta = sn@meta.data
}

# Check of name right
if(!opt$sntransfer %in% names(sn@meta.data)){
    warning(str_glue('{opt$sntransfer} not found in snRNA meta data. Please double check file'))
    warning(str_glue('Available ones are {toString(names(sn@meta.data))}'))
    stop("Please double check file")
}

# Clean name - prohibited character : "/"
sn@meta.data[[opt$sntransfer]] = str_replace(sn@meta.data[[opt$sntransfer]], '/','_')

# Run RCTD
RCTDresult = getRCTD_result(
    st_obj = ST,
    snrna_obj = sn, 
    colu_cell_type = opt$sntransfer,
    cell_count_filter = opt$cellthreshold
)

# Get sample name
if(is.null(opt$sample)){
    opt$sample = str_remove_all(opt$st, '.+/|[[.]]rds|[[.]]qs|[[.]]RDS')
}

# Process RCTD result
assay = RCTD2assay(RCTDresult, add_id = F, return_as = 'assay')
meta  = RCTD2meta(assay, col_name = str_glue('RCTD_{opt$sntransfer}'))
meta_newall=integrate_RCTDmeta(orig_meta =  ST@meta.data, new_meta = meta)

# ADD TO OBJ
ST@meta.data =  meta_newall
ST@assays[[str_glue('RCTD_{opt$sntransfer}')]] = assay

# Outputs:
# Save RCTD
dir.create(str_glue('{opt$out}/{opt$sample}'), recursive = T)
saveRDS(RCTDresult, str_glue('{opt$out}/{opt$sample}/RCTDresult_{opt$transfername}_{opt$sntransfer}.rds'))

# Save processed result
saveRDS(assay, str_glue('{opt$out}/{opt$sample}/RCTDassay_{opt$transfername}_{opt$sntransfer}.rds'))
saveRDS(meta_newall, str_glue('{opt$out}/{opt$sample}/RCTDmeta_{opt$transfername}_{opt$sntransfer}.rds'))

# SAVE final ST
message('Saving final ST obj')
saveRDS(ST, str_glue('{opt$out}/{opt$sample}/STobj_RCTD_{opt$transfername}_{opt$sntransfer}.rds'))
