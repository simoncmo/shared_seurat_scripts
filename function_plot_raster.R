# Rasterize layer
RasterizeGGLayer = function(p, raster_layers, dpi = 500){
    if(missing(raster_layers)) raster_layers = 1:length(p$layers)
    p$layers = imap(p$layers, function(layer, idx){
        if(idx %in% raster_layers) ggrastr::rasterise(layer, dpi = dpi) # Raster
        else layer # unchange
    })
    p
}