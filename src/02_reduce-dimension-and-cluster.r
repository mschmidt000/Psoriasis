### reduce dimension and find clusters
### 19.05.22

# reduce dimension and plot -----------------------------------------------

for(i in seq_along(runs)){
  seurat_object <- load_seurat_object(data_set_name = runs[i],  output_data_path = output_data_path) %>%
  reduce_dimension_and_cluster(figures_path = figures_path)
  
  save_seurat_object(seurat_object = seurat_object, output_data_path = output_data_path)
  rm(seurat_object)
}

