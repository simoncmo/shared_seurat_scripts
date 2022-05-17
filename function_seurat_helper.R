library(tidyverse)

# get color
get_rainbow_col = function(n=10, start=0.2, end =1  ){
col_rainbow = c("#E7EBFA", "#E7E9F9", "#E5E7F8", "#E4E4F6", "#E2E2F5", 
"#E1DFF3", "#DFDCF2", "#DEDAF0", "#DDD7EE", "#DBD4ED", "#D9D1EB", 
"#D8CEE9", "#D6CBE7", "#D5C8E5", "#D3C5E4", "#D2C2E2", "#D0BFE0", 
"#CEBCDE", "#CCB8DC", "#CAB5D9", "#C9B2D7", "#C7AFD5", "#C5ACD3", 
"#C3A8D1", "#C1A5CF", "#C0A2CD", "#BE9FCB", "#BC9BC9", "#BA98C7", 
"#B895C5", "#B692C4", "#B58EC2", "#B38BC0", "#B188BE", "#AF85BC", 
"#AD82BA", "#AC7FB9", "#AA7CB7", "#A879B5", "#A677B3", "#A574B1", 
"#A371B0", "#A26EAE", "#A06BAC", "#9E68AB", "#9D65A9", "#9B63A7", 
"#9960A6", "#975DA4", "#965BA2", "#9458A0", "#92559E", "#90539C", 
"#8E509B", "#8C4E99", "#884E99", "#844D99", "#814C99", "#7D4D99", 
"#794D99", "#754C9A", "#714C9A", "#6E4C9B", "#6D4D9D", "#6B4F9E", 
"#6951A0", "#6753A2", "#6554A4", "#6356A6", "#6158A7", "#5F59A9", 
"#5E5BAB", "#5D5DAD", "#5B5FAF", "#5A61B1", "#5863B3", "#5765B5", 
"#5567B7", "#5469B9", "#536BBA", "#536EBC", "#5270BE", "#5172BF", 
"#5074C1", "#4F77C3", "#4E79C4", "#4E7BC5", "#4E7DC5", "#4E7FC5", 
"#4E81C5", "#4E84C5", "#4D86C5", "#4D88C5", "#4D8AC5", "#4D8CC4", 
"#4D8DC3", "#4E8FC1", "#4E90C0", "#4E92BF", "#4E93BD", "#4E95BC", 
"#4E96BB", "#4F97BA", "#5098B9", "#5199B7", "#519AB6", "#529BB5", 
"#539CB4", "#539DB3", "#549EB1", "#559FB0", "#56A0AF", "#56A1AE", 
"#57A2AC", "#58A3AB", "#58A4AA", "#59A5A8", "#5AA5A7", "#5BA6A5", 
"#5CA7A4", "#5CA8A3", "#5DA8A1", "#5EA9A0", "#5FAA9E", "#60AB9C", 
"#61AC9B", "#63AC99", "#64AD97", "#65AE96", "#66AF94", "#67B092", 
"#68B091", "#69B18E", "#6BB28C", "#6DB289", "#6FB387", "#71B484", 
"#73B582", "#74B67F", "#76B77D", "#79B77A", "#7CB877", "#7FB875", 
"#82B972", "#84BA6F", "#87BA6D", "#89BB6A", "#8CBC67", "#90BC65", 
"#93BC62", "#97BD5F", "#9ABD5D", "#9EBD5A", "#A1BE58", "#A4BE55", 
"#A7BE53", "#AABE51", "#AEBD50", "#B1BD4E", "#B4BD4D", "#B7BC4B", 
"#BABC4A", "#BDBC48", "#BFBB47", "#C2BA46", "#C4B945", "#C7B844", 
"#C9B743", "#CCB642", "#CEB541", "#D0B440", "#D2B340", "#D4B23F", 
"#D5B03F", "#D7AF3E", "#D8AD3D", "#DAAC3D", "#DBAB3C", "#DDA93B", 
"#DEA73B", "#DFA53A", "#E0A43A", "#E1A239", "#E2A039", "#E29E39", 
"#E39C39", "#E49A38", "#E59838", "#E59637", "#E59437", "#E59236", 
"#E59036", "#E68E35", "#E68C34", "#E78A34", "#E68733", "#E78533", 
"#E68233", "#E68033", "#E57D32", "#E67B32", "#E57831", "#E67631", 
"#E57330", "#E5702F", "#E56D2F", "#E56A2E", "#E4672D", "#E4642D", 
"#E3612C", "#E35E2C", "#E25B2B", "#E1572A", "#E1542A", "#E05029", 
"#DF4D28", "#DF4927", "#DE4527", "#DD4126", "#DD3D25", "#DC3824", 
"#DC3424", "#DB2E23", "#DA2922", "#D92222", "#D52221", "#D12220", 
"#CD2220", "#C8221F", "#C4221F", "#BF221E", "#BB221E", "#B7221E", 
"#B2221D", "#AD221D", "#A9221C", "#A4221C", "#A0211C", "#9B211B", 
"#97211A", "#92211A", "#8E2019", "#892019", "#841F19", "#801F18", 
"#7B1F18", "#771E17", "#721E17", "#6E1D16", "#6A1D15", "#661C15", 
"#611C14", "#5D1B14", "#591A13", "#551A13", "#521913")

  idx = floor(seq(from = start, to = end, length.out = n)*256)
  return(col_rainbow[idx])
}

# Set spatial feature fill scale
spatial_feature_unify_fill_scale = function(plot_obj_list, st_obj,  markers){
  colors = colorRampPalette(c('#5953A3','#8ACFA4','#FEFDBC','#F57346','#A21743'))(100)
  nmarker = length(markers)
  limit_max =GetAssayData(st_obj)[markers, ] %>% as.array %>% apply(1, max)
  limit_min =GetAssayData(st_obj)[markers, ] %>% as.array %>% apply(1, min)
  plot_obj_list = imap(plot_obj_list, function(obj, idx){
    i = ifelse(idx %% nmarker!=0, idx %% nmarker, nmarker) 
    obj + scale_fill_gradientn(colors = colors, limits=c(limit_min[[i]],limit_max[[i]]))
  })
  plot_obj_list
}

# Rotate image 2/2/2021
# Helper s
 rotate_matrix_clockwise = function(mat, rotate=90){
   if(rotate==90){ # (x', y') = (y0, -x0)
     t(mat)[ ,ncol(t(mat)):1]
   }else if(rotate==270 || rotate ==-90){  # (x', y') = (-y0, x0)
     t(mat[,ncol(mat):1])
   }else if(rotate==180){ # (x', y') = (-x0, -y0)
     mat[nrow(mat):1, ncol(mat):1] 
   }else{
     stop("Only support rotate = 90, 180, -90(270) rotations.")
   }
 }

rotate_image = function (image, angle) 
{
  if (inherits(image, "data.frame")) {
    image = as.matrix(image)
  }else if(all(c(inherits(image, "array"), !is.na(dim(image)[3]), dim(image)[3] == 3))){
    new_array = map(1:dim(image)[3], function(i){
      rotate_matrix_clockwise(image[, , i], angle)
    })
    out = array(unlist(new_array), dim = c(nrow(new_array[[1]]),
                                           ncol(new_array[[1]]), length(new_array)))
  }else{
    stop("invalid type of image, supported types are matrix and 3 dimensional array")
  }
  return(out)
}

rotate_coord = function(coord,rotate, x_shift=0, y_shift=0){ # Flip x and y, then turn new x negative
  if(rotate==90){
    # Flip imagerow and imagecol value
    coord = setNames(coord, c('tissue', 'row','col','imagecol','imagerow')) 
    coord = coord[, c('tissue','row','col','imagerow','imagecol')] 
    #coord$col      = -coord$col + (max(coord$col)-min(coord$col)) # Add negative to new x(column)
    coord$imagecol = -coord$imagecol + (max(coord$imagecol)-min(coord$imagecol)) 
  }else if(rotate==270 || rotate ==-90){
    # Flip imagerow and imagecol value
    coord = setNames(coord, c('tissue', 'row','col','imagecol','imagerow')) 
    coord = coord[, c('tissue','row','col','imagerow','imagecol')] 
    #coord$row      = -coord$row + (max(coord$row)-min(coord$row)) # Add negative to new x(column)
    coord$imagerow = -coord$imagerow + (max(coord$imagerow)-min(coord$imagerow)) # Add negative to new x(column)
  }else if(rotate==180){ # add negative on both axis
    coord$imagecol = -coord$imagecol + (max(coord$imagecol)-min(coord$imagecol)) 
    coord$imagerow = -coord$imagerow + (max(coord$imagerow)-min(coord$imagerow)) 
  }else{
    stop("Only support rotate = 90, 180, -90(270) rotations.")
  }
  
  # shift the image
  coord$imagecol = coord$imagecol + x_shift
  coord$imagerow = coord$imagerow + y_shift
  print("Note, row and col in coordinates are not altered.")
  coord
}

# Main image rotation function
rotate_spatial_slice = function(obj, slice='slice1', rotate, x_shift=0, y_shift=0){
  if(missing(rotate) || !rotate %in% c(90,180,-90,270)){
    stop("Only support rotate = 90, 180, -90(270) rotation")
  }
  # Rotate coordinate
  coord = obj@images[[slice]]@coordinates
  obj@images[[slice]]@coordinates = rotate_coord(coord,rotate, x_shift, y_shift)
  
  # Rotate image
  obj@images[[slice]]@image = rotate_image(obj@images[[slice]]@image, rotate)
  print(str_glue("Rotated {rotate} degree. Adjust x_shift or y_shift if not aligned with dot."))
  obj
}

# Get 
getDominantAssignment=function(obj, assay_name){
  if(!assay_name%in%names(obj@assays)) {stop(str_glue("Assay {assay_name} not found."))}
  celltype_percent = obj@assays[[assay_name]] %>% GetAssayData %>% as.matrix %>% t %>% as.data.frame
  celltype_df     = data.frame(row.names = rownames(celltype_percent),
                               celltype  = colnames(celltype_percent)[max.col(celltype_percent, ties.method = 'first')])
  celltype_df
}


# Get marker scores from a list of gene
helper_getMarkerScore = function(obj, features){
  if(missing(features)) {stop("please provide features to average")}
  else{
    features = intersect(features, rownames(obj))
    message(str_glue("Using {length(features)} features."))
    GetAssayData(obj)[features, ] %>% apply(2, mean) 
  }
}

# Get marker scores from a gene dataframe
getAllMarkerScore = function(obj, marker_df){
  if(length(intersect(c('gene_symbol','cell_type'), colnames(marker_df)))!=2){
    stop("Column names for the marker df has to be gene_symbol and cell_tyoe")
  }else{
    lapply(unique(marker_df$cell_type), function(type){
      features = filter(marker_df, cell_type==type)$gene_symbol
      helper_getMarkerScore(obj, features)
    }) %>%
      bind_rows %>% 
      t %>% as.data.frame %>%
      setNames(unique(marker_df$cell_type))
  }
}

