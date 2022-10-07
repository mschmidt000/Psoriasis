### resort opossom run
### 25.08.22
pacman::p_load(oposSOM, here, Seurat)
library(oposSOM)
library(here)
library
source(here("src", "opossom-functions.r"))
source(here("src", "paths.r"))

filename <- here(output_data_path, "integrated-seurat-obj.RData")
load(filename)
filename <- here("2022-08-24_psoriasis.RData")
load(filename)

IL18_data <- GetAssayData(obj_integr, slot = "data")["IL18",]
obj_integr$IL18up <- FALSE
obj_integr$IL18up[names(IL18_data[which(IL18_data > (mean(IL18_data) + sd(IL18_data)))])] <- TRUE

IL32_data <- GetAssayData(obj_integr, slot = "data")["IL32",]
obj_integr$IL32up <- FALSE
obj_integr$IL32up[names(IL32_data[which(IL32_data > (mean(IL32_data) + sd(IL32_data)))])] <- TRUE

obj_integr$predicted.id.Cell_type_nicknames_IL18_IL32 <- as.character(obj_integr$predicted.id.Cell_type_nicknames)
obj_integr$predicted.id.Cell_type_nicknames_IL18_IL32[which(obj_integr$IL18up)] <- 
  paste(obj_integr$predicted.id.Cell_type_nicknames_IL18_IL32[which(obj_integr$IL18up)], "IL18up", sep = "_")
obj_integr$predicted.id.Cell_type_nicknames_IL18_IL32[which(obj_integr$IL32up)] <- 
  paste(obj_integr$predicted.id.Cell_type_nicknames_IL18_IL32[which(obj_integr$IL32up)], "IL32up", sep = "_")

env <- resort_env(env, obj_integr$predicted.id.Cell_type_nicknames_IL18_IL32, "IL18up-IL32up-subpopulations")
analysis_and_output_after_resort(env)
