### run exSOM
### 2022-08-24_psoriasis_predicted-celltype + billi2022_seurat-object
### 07.09.22


# load packages and source code -------------------------------------------

pacman::p_load(Seurat, tidyverse, oposSOM, here)
source(here("src","seurat-functions.r"))
source(here("src","opossom-functions.r"))
source(here("src","paths.r"))

# load data ---------------------------------------------------------------

filename <- here("2022-08-24_psoriasis_predicted-celltype.RData")
load(filename)

filename <- here(output_data_path, "billi2022_seurat-object.RData")
load(filename)
Idents(seurat_object) <- predicted.id.Cell_type_nicknames
seurat_object <- subset(seurat_object, downsample = 100)




# run programm ------------------------------------------------------------

env <- run_exSOM(env, seurat_object, "reynolds2021" )
