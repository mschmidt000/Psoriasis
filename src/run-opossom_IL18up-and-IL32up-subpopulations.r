### run opossom on melanocyte cells
### 26.05.22
library(oposSOM)

filename <- here(output_data_path, "integrated-seurat-obj.RData")
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

Idents(obj_integr) <- "predicted.id.Cell_type_nicknames_IL18_IL32"
# obj_mel_ds <- subset(obj_mel, downsample = 300) %>%
#   ScaleData(verbose = FALSE, features = rownames(obj_mel))
# rm(obj_mel)
env <- opossom.new(list(
  dataset.name = paste(Sys.Date(), "psoriasis", sep = "_"),
  dim.1stLvlSom = "auto",
  dim.2ndLvlSom = "auto",
  activated.modules = list(
    "reporting" = TRUE,
    "primary.analysis" = TRUE,
    "sample.similarity.analysis" = FALSE,
    "geneset.analysis" = TRUE,
    "psf.analysis" = TRUE,
    "group.analysis" = TRUE,
    "difference.analysis" = TRUE
  ),
  standard.spot.modules = "group.overexpression"
))


# definition of indata, group.labels and group.colors
env$indata <- GetAssayData(obj_integr, slot = "scale.data", assay = "integrated") %>%
  as.matrix()
env$group.labels <- Idents(obj_integr)
env$group.colors <- env$group.labels %>%
  n_distinct() %>%
  scales::hue_pal()(.)
names(env$group.colors) <- levels(env$group.labels) 
env$group.colors <- env$group.colors[match(env$group.labels, names(env$group.colors))]
names(env$group.colors) <- names(env$group.labels)
ord <- order(env$group.labels)
env$indata <- env$indata[, ord]
env$group.labels <- env$group.labels[ord]
env$group.colors <- env$group.colors[ord]
opossom.run(env)

