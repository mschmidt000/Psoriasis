### receptor-ligand analysis for every single data set and the complete integrated data set
### 25.05.22

library(liana)
source(here("src", "liana-functions.r"))
filename <- here(output_data_path, "integrated-seurat-obj.RData")
load(filename)


Idents(obj_integr) <- "orig.ident"
obj_naevus_neg <- subset(obj_integr, idents = "LE497NA_Rep", invert = TRUE)
Idents(obj_naevus_neg) <- interesting_idents_clust_mels

liana_obj <- liana_wrap(obj_integr, assay = "RNA")
filename <- here(output_data_path, "integrated-liana-obj.RData")
save(liana_obj, file = filename)
liana_obj <- liana_obj %>%
  liana_aggregate()
save(liana_obj, file = filename)

plot_lymphs_vs_mels(liana_obj, obj_naevus_neg, "integrated")

Idents(obj_integr) <- "orig.ident"

sapply(levels(Idents(obj_integr)), function(x) {
  obj_subset <- subset(obj_integr, idents = x)
  obj_subset <- ScaleData(obj_subset, verbose = FALSE, features = rownames(obj_subset)) %>%
    RunPCA(npcs = n_dims_use, verbose = FALSE)

  liana_obj <- liana_wrap(obj_subset, assay = "RNA")
  filename <- here(output_data_path, paste0(x, "-liana-obj.RData"))
  save(liana_obj, file = filename)
  liana_obj <- liana_obj %>%
    liana_aggregate()
  save(liana_obj, file = filename)

  plot_lymphs_vs_mels(liana_obj, obj_subset, x)
})
