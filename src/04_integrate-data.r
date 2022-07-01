### integrating good quality data sets
### 20.05.22
source(here("src", "make-nicknames-for-celltypes.r"))

n_dims_use <- 30

runs_to_integrate <- c(
  "BF-LE-01-KT-PSO_all",  "BF-LE-02-PG-PSO_all",  "BF_LE_03_VC_03_al"
) %>% sort()

runs_list <- map(runs_to_integrate, ~load_seurat_object(data_set_name = .x,  output_data_path = output_data_path))
names(runs_list) <- runs_to_integrate
integr_features <- map( runs_list, ~rownames(.))
integr_features <- Reduce(intersect, integr_features)
integr_anchors <- FindIntegrationAnchors(runs_list, dims = 1:n_dims_use, anchor.features = integr_features)

rm(runs_list)

obj_integr <- IntegrateData(anchorset = integr_anchors, dims = 1:30)
DefaultAssay(obj_integr) <- "integrated"
filename <- here(output_data_path, "integrated-seurat-obj.RData")
save(obj_integr, file = filename)
rm(integr_anchors)

obj_integr <- ScaleData(obj_integr, verbose = FALSE, features = rownames(obj_integr)) %>%
  RunPCA(npcs = n_dims_use, verbose = FALSE) %>%
  RunTSNE(reduction = "pca", dims = 1:n_dims_use, perplexity = sqrt(ncol(obj_integr))) %>%
  RunUMAP(dims = 1:n_dims_use) %>%
  FindNeighbors(reduction = "pca", dims = 1:n_dims_use) %>%
  FindClusters(resolution = c(0.2, 0.4, 0.5, 0.6, 0.8, 1, 1.2))

obj_integr$predicted.id.Cell_type_nicknames <- make_nicknames(obj_integr$predicted.id.Cell_type)
obj_integr$predicted.id.Cell_type <- factor(obj_integr$predicted.id.Cell_type, levels = sort(unique(obj_integr$predicted.id.Cell_type)))
obj_integr$predicted.id.Cell_group <- factor(obj_integr$predicted.id.Cell_group, levels = sort(unique(obj_integr$predicted.id.Cell_group)))
obj_integr$predicted.id.Flow_gate <- factor(obj_integr$predicted.id.Flow_gate, levels = sort(unique(obj_integr$predicted.id.Flow_gate)))

filename <- here(output_data_path, "integrated-seurat-obj.RData")
save(obj_integr, file = filename)

interesting_idents_long <- c("predicted.id.Cell_type_nicknames", "predicted.id.Cell_group", "predicted.id.Flow_gate")
interesting_idents_short <- c("Reynolds 21 Cell Type", "Reynolds 21 Cell Group", "Reynolds 21 Flow Gate")

color_list <- list()
for (i in interesting_idents_long) {
  color_list[[i]] <- randomcoloR::distinctColorPalette(n_distinct(obj_integr@meta.data[, i]))
  names(color_list[[i]]) <- levels(obj_integr@meta.data[, i])
}
filename <- here(output_data_path, "color-list.RData")
save(color_list, file = filename)

filename <- here(figures_path, "04_integrated-dimension-reduction-and-clustering.pdf")
pdf(filename, height = 7, width = 14)
DimPlot(obj_integr, reduction = "umap", label = TRUE, group.by = "orig.ident", repel = TRUE) + ggtitle("UMAP colored by Run")
DimPlot(obj_integr, reduction = "umap", label = FALSE, group.by = "orig.ident", repel = TRUE) + ggtitle("UMAP colored by Run")
barplot(table(obj_integr@meta.data$orig.ident), main = paste0("Number of Cells per Sample"), col = scales::hue_pal()(length(table(obj_integr@meta.data$orig.ident))))
my_plot <- DimPlot(obj_integr, reduction = "umap", label = TRUE, group.by = paste0("integrated_snn_res.", c(0.5, 0.6, 0.8, 1)))
my_plot
clustree::clustree(x = obj_integr@meta.data, prefix = "integrated_snn_res.")
DimPlot(obj_integr, reduction = "umap", label = FALSE, repel = TRUE, group.by = "Phase")
my_plot <- FeaturePlot(obj_integr, reduction = "umap", order = TRUE, features = c("percent.mt", "percent.rps", "nCount_RNA", "nFeature_RNA"))
my_plot
DotPlot(obj_integr, features = c("percent.mt", "percent.rps"), col.min = 0, col.max = 100, dot.scale = 4)
dev.off()


Idents(obj_integr) <- "orig.ident"
filename <- here(figures_path, "04_integrated-dimension-reduction-col-by-predicted-cell-types.pdf")
pdf(filename, height = 7, width = 14)
plot(0, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(0, 1))
mtext("All 11 Data Sets", side = 3, line = -2, cex = 3, at = -0.04, font = 3, adj = 0)

for (i in seq_along(interesting_idents_long)) {
  obj_integr@meta.data[, interesting_idents_long[i]] <- factor(obj_integr@meta.data[, interesting_idents_long[i]], level = sort(unique(obj_integr@meta.data[, interesting_idents_long[i]])))
  my_title <- paste0("UMAP colored by ", interesting_idents_short[i])
  print(DimPlot(obj_integr, reduction = "umap", label = TRUE, group.by = interesting_idents_long[i], repel = TRUE, cols = color_list[[i]]) + ggtitle(my_title))
  print(DimPlot(obj_integr, reduction = "umap", label = FALSE, group.by = interesting_idents_long[i], repel = TRUE, cols = color_list[[i]]) + ggtitle(my_title))
  par(mar = c(10, 5, 4, 5))
  my_title <- paste0("Number of Cells per Cluster (", interesting_idents_short[i], ")")
  my_table <- table(obj_integr@meta.data[, interesting_idents_long[i]])
  print(barplot(my_table, las = 2, main = my_title, col = color_list[[i]][names(my_table)]))
  rm(my_table)
}

for (a in seq_along(levels(Idents(obj_integr)))) {
  plot(0, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(0, 1))
  mtext(levels(Idents(obj_integr))[a], side = 3, line = -2, cex = 3, at = -0.04, font = 3, adj = 0)
  seurat_obj <- subset(obj_integr, idents = levels(Idents(obj_integr))[a])
  for (i in seq_along(interesting_idents_long)) {
    my_title <- paste0("UMAP colored by ", interesting_idents_short[i], " (", levels(Idents(obj_integr))[a], ")")
    print(DimPlot(seurat_obj, reduction = "umap", label = TRUE, group.by = interesting_idents_long[i], repel = TRUE, cols = color_list[[i]][unique(seurat_obj@meta.data[, interesting_idents_long[i]])]) + ggtitle(my_title))
    print(DimPlot(seurat_obj, reduction = "umap", label = FALSE, group.by = interesting_idents_long[i], repel = TRUE, cols = color_list[[i]][unique(seurat_obj@meta.data[, interesting_idents_long[i]])]) + ggtitle(my_title))
    par(mar = c(10, 5, 4, 5))
    my_title <- paste0("Number of Cells per Cluster (", interesting_idents_short[i], ", ", levels(Idents(obj_integr))[a], ")")
    my_table <- table(seurat_obj@meta.data[, interesting_idents_long[i]])
    print(barplot(my_table, las = 2, main = my_title, col = color_list[[i]][names(my_table)]))
  }
}
dev.off()
