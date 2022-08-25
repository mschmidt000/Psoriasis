### run opossom on melanocyte cells
### 26.05.22
library(oposSOM)

filename <- here(output_data_path, "integrated-seurat-obj.RData")
load(filename)

annotated_clusters <- annotate_clusters(obj_integr, "predicted.id.Cell_type_nicknames")
obj_integr$seurat_clusters_annotated <- annotated_clusters[obj_integr$seurat_clusters]
obj_integr$seurat_clusters_annotated <- paste0(obj_integr$seurat_clusters," (",obj_integr$seurat_clusters_annotated , ")")
obj_integr$seurat_clusters_annotated <- factor(obj_integr$seurat_clusters_annotated, levels = unique(obj_integr$seurat_clusters_annotated)[order(unique(obj_integr$seurat_cluster))])

Idents(obj_integr) <- "seurat_clusters_annotated"
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

