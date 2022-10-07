### resort opossom run
### 25.08.22
library(oposSOM)
library(here)
library(Seurat)
source(here("src", "opossom-functions.r"))
source(here("src", "paths.r"))
filename <- here(output_data_path, "color-list.RData")
load(filename)  


filename <- here(output_data_path, "integrated-seurat-obj.RData")
load(filename)
filename <- here("2022-08-24_psoriasis.RData")
load(filename)

new_labels <- as.character(obj_integr$predicted.id.Cell_type_nicknames)
names(new_labels) <- rownames(obj_integr@meta.data)

env <- resort_env(env, new_labels, "predicted-celltype", color_list[["predicted.id.Cell_type_nicknames"]])
analysis_and_output_after_resort(env)

