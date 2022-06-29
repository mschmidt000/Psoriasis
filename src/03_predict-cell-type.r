### predict cell types with data set from reynolds 2021
### 20.05.22

cell_type_to_exclude <- "Schwann"

filename <- here(literature_path, "Reynolds2021_Data-all_seuratObj.RData")
load(filename)
Idents(seurObj) <- "Cell_type"
# seurObj <- subset(seurObj, idents = cell_type_to_exclude, invert = TRUE)

for(i in seq_along(runs)){
  seurat_object <- load_seurat_object(data_set_name = runs[i],  output_data_path = output_data_path) %>%
                       transfer_labels(reference_object = reference_object, figures_path = figures_path)
  
  save_seurat_object(seurat_object = seurat_object, output_data_path = output_data_path)
  rm(seurat_object)
}
