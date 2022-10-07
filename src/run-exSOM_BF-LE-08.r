### run exSOM
### 2022-08-24_psoriasis_predicted-celltype + BF_LE_08_GD_pre_h
### 06.09.22


# load packages and source code -------------------------------------------

pacman::p_load(Seurat, tidyverse, oposSOM, here)
source(here("src","seurat-functions.r"))
source(here("src","opossom-functions.r"))
source(here("src","paths.r"))

# load data ---------------------------------------------------------------

filename <- here("2022-08-24_psoriasis_predicted-celltype.RData")
load(filename)

seurat_object <- load_seurat_object(data_set_name = "BF_LE_08_GD_pre_h",  output_data_path = output_data_path)
Idents(seurat_object) <- "predicted.id.Cell_type_nicknames"
seurat_object <- subset(seurat_object, downsample = 400)
  
# run programm ------------------------------------------------------------

env <- run_exSOM(env, seurat_object, "BF_LE_08_GD_pre_h")

