
# Wrapper to load Viisum from Spaceranger output
Load10X_STutility = function(data_dir, ...){
  # Load STutility object
  # data_dir: path to the 10x output folder
  # slice: slice name
  # return: STutility object
  # choose tissue position file name 
  if(file.exists(str_glue('{data_dir}/tissue_positions_list.csv'))){
    tissue_positions_file = 'tissue_positions_list.csv'
    }else{
    tissue_positions_file = 'tissue_positions.csv'
    }
  infotable = data.frame(samples = file.path(data_dir, 'filtered_feature_bc_matrix.h5'),
    spotfiles = file.path(data_dir, 'spatial', tissue_positions_file),
    imgs = file.path(data_dir, 'spatial', 'tissue_hires_image.png'),
    json = file.path(data_dir, 'spatial', 'scalefactors_json.json'))
  
  # Load STutility object
  STutility::InputFromTable(infotable, platform = 'Visium', ...)
}