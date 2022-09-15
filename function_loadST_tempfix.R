############# ADHOC-fix ########################
################################################################
message('Temprary fix for Seurat Load10X_Spatial function till they release update')
message('Use Load10X_Spatial_fixed to load the ST file')
# Loading a bit more complicated now. Need fix. or update Seurat
Load10X_Spatial_fixed = function (data.dir, filename = "filtered_feature_bc_matrix.h5", 
          assay = "Spatial", slice = "slice1", filter.matrix = TRUE, 
          to.upper = FALSE, ...) 
{
  if (length(x = data.dir) > 1) {
    warning("'Load10X_Spatial' accepts only one 'data.dir'", 
            immediate. = TRUE)
    data.dir <- data.dir[1]
  }
  data <- Read10X_h5(filename = file.path(data.dir, filename), 
                     ...)
  if (to.upper) {
    rownames(x = data) <- toupper(x = rownames(x = data))
  }
  object <- CreateSeuratObject(counts = data, assay = assay)
  image <- Read10X_Image_fixed(image.dir = file.path(data.dir, "spatial"), 
                         filter.matrix = filter.matrix)
  image <- image[Cells(x = object)]
  DefaultAssay(object = image) <- assay
  object[[slice]] <- image
  return(object)
}

Read10X_Image_fixed= function (image.dir, image.name = "tissue_lowres_image.png", 
                                  filter.matrix = TRUE, ...) 
{
  image <- png::readPNG(source = file.path(image.dir, image.name))
  scale.factors <- jsonlite::fromJSON(txt = file.path(image.dir, "scalefactors_json.json"))
  tissue.positions <- read.csv(file = file.path(image.dir, 
                                                "tissue_positions.csv"), col.names = c("barcodes", 
                                                                                            "tissue", "row", "col", "imagerow", "imagecol"), header = TRUE, 
                               as.is = TRUE, row.names = 1)
  if (filter.matrix) {
    tissue.positions <- tissue.positions[which(x = tissue.positions$tissue == 
                                                 1), , drop = FALSE]
  }
  unnormalized.radius <- scale.factors$fiducial_diameter_fullres * 
    scale.factors$tissue_lowres_scalef
  spot.radius <- unnormalized.radius/max(dim(x = image))
  return(new(Class = "VisiumV1", image = image, scale.factors = scalefactors(spot = scale.factors$tissue_hires_scalef, 
                                                                             fiducial = scale.factors$fiducial_diameter_fullres, hires = scale.factors$tissue_hires_scalef, 
                                                                             scale.factors$tissue_lowres_scalef), coordinates = tissue.positions, 
             spot.radius = spot.radius))
}
