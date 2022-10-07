### preprocess the reference data set of lupus erythematodes cells
### 02.09.22

# load packages and source code -------------------------------------------

pacman::p_load(Seurat, tidyverse, clustree)
source(here("src","seurat-functions.r"))
source(here("src","paths.r"))
source(here("src", "make-nicknames-for-celltypes.r"))


# set variables and load data ---------------------------------------------

input_data_path <- here("literature", "GSE186476")
paths <- runs <- list.files(input_data_path)
n_dims_use <- 30

filename <- here(literature_path, "Reynolds2021_Data-all_seuratObj.RData")
load(filename)
reference_object <- seurat_object

# preprocess, filter and plot ---------------------------------------------
runs_list <- list()

for(i in 4:length(runs)){
  seurat_object <- create_seurat_object(  
    data_set_name = runs[i],
    input_data_folder = paths[i],
    input_data_path = input_data_path
  ) %>%
    NormalizeData() %>%
    FindVariableFeatures(do.plot = FALSE, verbose = FALSE) %>%
    ScaleData() %>%
    reduce_dimension_and_cluster(figures_path = figures_path) %>%
    transfer_labels(reference_object = reference_object, figures_path = figures_path)
  
  save_seurat_object(seurat_object, output_data_path = output_data_path)
  runs_list[[runs[i]]] <- seurat_object
  rm(seurat_object)
}

integr_features <- map( runs_list, ~rownames(.))
integr_features <- Reduce(intersect, integr_features)
integr_anchors <- FindIntegrationAnchors(runs_list, dims = 1:30, anchor.features = integr_features)

rm(runs_list)

seurat_object <- IntegrateData(anchorset = integr_anchors, dims = 1:30)
DefaultAssay(seurat_object) <- "integrated"
filename <- here(literature_path, "billi2022_seurat-object.RData")
save(seurat_object, file = filename)
rm(integr_anchors)

seurat_object <- ScaleData(seurat_object, verbose = FALSE, features = rownames(seurat_object)) %>%
  RunPCA(npcs = n_dims_use, verbose = FALSE) %>%
  RunTSNE(reduction = "pca", dims = 1:n_dims_use, perplexity = sqrt(ncol(seurat_object))) %>%
  RunUMAP(dims = 1:n_dims_use) %>%
  FindNeighbors(reduction = "pca", dims = 1:n_dims_use) %>%
  FindClusters(resolution = c(0.2, 0.4, 0.5, 0.6, 0.8, 1, 1.2))


filename <- here(literature_path, "billi2022_seurat-object.RData")
save(seurat_object, file = filename)
filename <- here(output_data_path, "billi2022_seurat-object.RData")
save(seurat_object, file = filename)

