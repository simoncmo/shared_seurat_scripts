

################################################
# Make data set
################################################
GetKidenyDataSet = function(dataset = c('Human_healthy_adult_kidney','Human_healthy_adult_kideny_epithilia')){
    dataset = match.arg(dataset)

    ################################################
    # Human CellType marker
    # https://humphreyslab.com/SingleCell/displaycharts.php
    # 1. Healthy Adult Kidney - Epithilia: Wu and Malone et al, JASN 2018.
    # Epi only
    cell_type_markers_epi = 
    list(CD_IC=c("ADGRF5","SLC4A1","KIT","AC015910.1","CLNK","RHCG","SLC26A7","NXPH2","PDE1C","IGFBP5","AC112229.3","SPINK1","AP000322.1","TMEM213","FZD7","AC026333.3","STAP1","RCAN2","SEMA3C"),
    CD_PC=c("AQP2","STC1","BMPR1B","HSD11B2","GATA3","PIK3C2G","PWRN1","MECOM","COBLL1","SLC38A1","AQP3","ADAMTS16","MPPED2","CADPS2","CRYBG1","L1CAM","TOX","TACSTD2","LHX1","TEX41"),
    DT=c("SLC12A3","TRPM6","AC019197.1","SLC8A1","TEX41","CNNM2","KLHL3","SCN2A","EFNA5","KCNIP4","CACNB4","ERBB4","ESRRG","EGF","KITLG","FMN1","PIK3C2G","CPEB4","TRPM7","MECOM"),
    LH=c("SLC12A1","ERBB4","GPC5","UMOD","PLCB1","PDE1A","PCDH9","KCNIP4","LINC01762","MECOM","RP1","ACPP","ESRRG","CASR","ENOX1","EGF","ESRRB","CCSER1","PKP4","HS6ST2"),
    P=c("PTPRQ","PTPRO","CLIC5","NPHS2","AC092813.2","ST6GALNAC3","ATP10A","NTNG1","PLA2R1","FYN","CTGF","LINC02149","SRGAP2C","PLCE1","ALS2CL","KIRREL1","PODXL","CDC42EP3","CDC14A","SRGAP2B"),
    PT=c("NEAT1","LRP2","NRP1","LINC00621","ACSM2A","BNC2","ACSM2B","C1orf186","SLC13A3","LINC01320","NLGN1","SNX29","MAST4","SYNE2","SASH1","WDR72","PTPRD","RAB11FIP3","NCKAP5","CUBN"))
    ################################################

    ################################################
    # Human CellType marker
    # https://humphreyslab.com/SingleCell/displaycharts.php
    # 3. Healthy Adult Kidney - Complete: Wu and Uchimura et al, Cell Stem Cell 2018; and RBK RID: 14-4KPM.
    cell_type_markers = list(Podocyte=c("NTNG1","CLIC5","NPHS1","FMN2","NPHS2","ADAMTS19","TARID","ATP10A","AC092813.2","CR1","ZNF804A","WT1","ALS2CL","NFASC","PCOLCE2","ARHGEF26","LINC02149","FGF1","WDR49","AC244021.1"),
    Mesangium=c("CARMN","EBF1","PDGFRB","NTRK3","ITGA8","SYNPO2","ADIRF","LMOD1","MYOCD","MICAL2","DAAM2","ADCY3","MEIS2","COL25A1","DGKG","CPED1","KALRN","TMTC1","RBMS3","NOTCH3"),
    EC=c("MEIS2","PTPRB","FLT1","NOTCH4","PLAT","EGFL7","AC010737.1","ERG","SLCO2A1","RAPGEF4","EMCN","BTNL9","PECAM1","PREX2","RUNX1T1","CDH13","KLF2","GFOD1","LDB2","ITGA8"),
    Macrophage=c("CD69","PIK3R5","THEMIS","ITK","CD247","CARD11","PTPRC","ARHGAP15","GPR141","BCL11B","PRKCB","DEF6","ADGRE5","STAT4","CHST11","CELF2","FLI1","INPP5D","POU2F2","AOAH"),
    'PT(S1)'=c("AL031599.1","AC105094.2","RNF212B","PCDH15","CDH20","ACSM2B","ACSM2A","SLC36A2","SLC5A12","SLC16A9","CYP3A5","LINC01060","CUBN","NPL","AP003400.1","NOX4","AFM","PAH","SLC4A4","PPP2R2B"),
    'PT(S2)'=c("ACSM2A","CUBN","ACSM2B","LRP2","AC096577.1","SLC13A3","SORCS1","SLC5A12","SLC4A4","TINAG","SLC16A9","SLC34A1","RAB11FIP3","HNF4A","AFM","SLC22A6","AC004053.1","SLC17A1","AC087762.1","ZNF804B"),
    'PT(S3)'=c("C5orf67","AL078590.3","SLC7A13","SLC30A8","ANKS1B","AC084128.1","GPM6A","SLC22A24","AL049629.1","LINC01789","ABCC3","AC128707.1","AL078590.2","SLC22A7","MOGAT1","CDH9","SYNE2","ZNF804B","C10orf126","AC098829.1"),
    'LH(DL)'=c("CFH","LINC01435","SLC4A11","LINC00924","KIRREL3","KLRG2","CLDN1","LINC01197","MSC-AS1","ALDH1A2","KCNT2","AC098617.1","GRM5","MGAT4C","FBXL7","HIST1H2AC","TSHZ2","ZFPM2","AC003984.1","SLIT3"),
    'LH(AL_1)'=c("AC092078.2","SLC12A1","CABP1","UMOD","PLCB1","CACNA2D3","GP2","SGIP1","ACPP","CALCR","CXCL12","CASR","ERBB4","SIM2","SCHLAP1","LINC01606","PHACTR1","EGF","ESRRB","GPC5"),
    'LH(AL_2)'=c("TMEM207","CLDN16","CASR","PLCB1","SLC12A1","ENOX1","UMOD","GP2","PRKD1","HIP1","ERBB4","LINC01762","SIM2","CACNA2D3","NPY1R","AC068631.1","FGF14","ARHGAP6","INPP4B","CCSER1"),
    DCT=c("SLC12A3","AC078980.1","TRPM6","KLHL3","LINC01055","TMEM52B","DEPDC1B","FMN1","CACNB4","SALL3","PID1","TBC1D9","CNNM2","ITPKB","TCF24","WNK1","TRIM50","EGF","TRPM7","LHX1"),
    CNT=c("LINC01098","SLC8A1","SLC8A1-AS1","TEX41","PDE4D","SCNN1G","SLC7A1","PWRN1","CADPS2","ABTB2","PRKG1","SNTG1","KCNIP1","ATP1B3","SLC38A11","AC019197.1","HSD11B2","AC015522.1","MYO1B","ARL15"),
    PC=c("SCNN1G","FXYD4","PWRN1","ROR2","SCNN1B","FGF12","CHGB","ZNF331","SLC7A1","PKIA-AS1","SLC9A4","EYA4","HUNK","NR4A2","GATA3","RASAL1","LINC01482","BMPR1B","ATP1B3","SLC38A1"),
    IC_A=c("SLC26A7","CLNK","ADGRF5","AC112229.1","HS6ST3","LINC00116","NXPH2","LY86-AS1","AC112229.3","AC092422.1","KIT","AC026333.3","STAP1","LEF1","HCAR2","DMRT2","AQP6","GALNT17","TMEM101","AC084048.1"),
    IC_B=c("SLC26A4","SAMHD1","SLC4A9","CA8","SLC35F3","NTRK1","TLDC2","AL136369.2","AC090502.1","SFTA1P","AC008438.1","CELF2","DGKI","PPM1E","AL355472.3","CLNK","AC007563.2","ATP6V1C2","PDE1C","LINC01187"),
    Undefined_1=c("RTP2","PAX3","PGM5P3-AS1","UMODL1","WNT3","AVPI1","DONSON","GABRE","USP17L7","LINC00290","HRNR","LINC01170","LINC01527","MT-ND5","MED12L","CCDC94","SMOC1","GRIA3","FP236383.2","LINC01340"),
    Undefined_2=c("LINC00621","TMEM178B","C1orf186","DCC","MEG3","PRICKLE1","MSC-AS1","CLSTN2","APBB1IP","CREB5","TPM1","ITGB8","ALPK2","MYO3A","UGT2B11","EGOT","VCAM1","COL23A1","PLA2G4C","EPB41L2"))
    ################################################

    # Reduce
    markers_more_than_one_cell_types = cell_type_markers %>% unlist %>% table %>% .[.>1] %>% names
    cell_type_markers_unique = map(cell_type_markers, ~setdiff(., markers_more_than_one_cell_types)) %>% .[!is.na(.)]


    if(dataset == 'Human_healthy_adult_kidney'){
        marker_use = cell_type_markers
    }else{
        marker_use = cell_type_markers_epi
    }

    # Reduce
    markers_more_than_one_cell_types = marker_use %>% unlist %>% table %>% .[.>1] %>% names
    cell_type_markers_unique = map(marker_use, ~setdiff(., markers_more_than_one_cell_types)) %>% .[!is.na(.)]

    # Plot
}


MakeKidneyMarkerDotPlot = function(obj, group.by = 'seurat_clusters', dataset = 'Human_healthy_adult_kidney', ...){
    marker_use = GetKidenyDataSet(dataset)
    DotPlot(obj, group.by = group.by, features = marker_use)+ RotatedAxis() + 
    labs(title = dataset) + 
    scale_color_gradientn(colors = RColorBrewer::brewer.pal(9, 'YlOrRd')[2:9]) + 
    theme(axis.text.x = element_text(size = 6))
}
