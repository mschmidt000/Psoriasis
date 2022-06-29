#' Make nicknames for predicted Reynolds et al. 2021 cell-types
#'
#' @param cell_types long Reynolds 2021 labels
#' @return short labels
make_nicknames <- function(cell_types) {
  cell_types <- gsub("Keratinocyte", "KC", cell_types)
  cell_types <- gsub("Macrophage", "Mac", cell_types)
  cell_types <- gsub("Fibroblast", "Fb", cell_types)
  cell_types <- gsub("Vascular endothelium", "VE", cell_types)
  cell_types <- gsub("Lymphatic endothelium", "LE", cell_types)
  cell_types <- gsub("Vascular endothelium", "VE", cell_types)
  cell_types <- gsub("mitotic", "mit", cell_types)
  cell_types <- gsub("Migratory", "Mig.", cell_types)
  cell_types <- factor(x = cell_types, levels = c(
       "KC premit", "KC mit", "KC postmit", "KC CD83+", "Melanocyte", "Stromal Schwann",
       "Fb 1", "Fb 2", "Fb 3", "Pericyte", "VE 1", "VE 2", "VE 3", "LE 1", "LE 2", "ILC", "NK1", "NK2", "NK3",
       "Tc", "Th", "Treg", "Mast cell", "Plasma cell", "Mac 1", "Mac 2", "DC1", "DC2", "LC1", "LC2","LC3", "LC4", "KLF10 LC",
       "Monocyte", "IL23 DC", "pDC", "moDC1", "moDC2", "Mig. cDC"
  ))
  cell_types
}



