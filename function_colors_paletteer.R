## Paletteer helper
library(paletteer)

# Function parameter = column name : https://stackoverflow.com/questions/62082563/function-to-filter-tibble-with-argument-having-same-name-as-a-column
ShowDiscretePalettes = function(package = 'Redmonder'){
	if(!package %in% palettes_d_names$package) stop("package not found. Please use palettes_d_names function to check available ones")

	# get all paletter under a package 
	palettes_d_names %>% filter(package == !!package) %>% 
		{setNames(pmap(., function(package, palette, ...) str_glue('{package}::{palette}') ), .$palette)} %>% map(paletteer_d)

} 

ExtendColorVectorPaletteer = function(color_vector = NULL, target_class, reduce = T, colors_use=NULL, randomize = F){
    missing_classes = setdiff(unique(target_class), names(color_vector))
	# Get color
	new_col_all = if(is.null(colors_use)){
		c(paletteer_d('Redmonder::qMSOBuWarm'),
		  #paletteer_d("Redmonder::qMSO12")[-2:0], # remove first 2 colors
		  #paletteer_d("Redmonder::qMSO15")[-2:0], # remove first colors
		  paletteer_d("Redmonder::qMSOGn")[-2:0], # remove first 2 colors
		  paletteer_d("Redmonder::qMSORdOr")[-2:0], # remove first 2 colors
		  paletteer_d("Redmonder::qMSORdPu")[-2:0] # remove first 2 colors

		  )
		}else{ colors_use }
    
    # rm duplicated color
    new_col_use = setdiff(new_col_all, color_vector) # remove exisiting color
    
    # Choose color
    new_col_use = if(randomize) sample(new_col_use,length(missing_classes)) else new_col_use[1:length(missing_classes)]
    new_col = new_col_use %>% setNames(missing_classes)
    
    # Combine
    color_vector = c(color_vector, new_col)
    
    # Remove 'Non exisiting item'toc)
    if(reduce) color_vector = color_vector[unique(target_class)]
}
