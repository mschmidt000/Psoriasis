### load data and preprocess with seurat
### 19.05.22


# load packages and source code -------------------------------------------

pacman::p_load(Seurat, tidyverse)
source(here("src","seurat-functions.r"))


# define variables --------------------------------------------------------

runs <- c(
  "BF-LE-01-KT-PSO_all",  "BF-LE-02-PG-PSO_all",  "BF_LE_03_VC_03_all",  "BF_LE_06_KS_LE_all",  "BF_LE_08_GD_pre_h"
)
input_data_path <- file.path("data")
output_data_path <- file.path("output")
figures_path <- file.path("figs")
literature_path <- "literature"

paths <- list.files(input_data_path, full.names = FALSE)


# preprocess, filter and plot ---------------------------------------------
for(i in seq_along(runs)){
  seurat_object <- create_seurat_object(  
                      data_set_name = runs[i],
                      input_data_folder = paths[i],
                      input_data_path = input_data_path
                ) %>%
                do_qc(figures_path = figures_path)  %>%
                do_filtering_and_qc(figures_path = figures_path)
              
  save_seurat_object(seurat_object, output_data_path = output_data_path)
  rm(seurat_object)
}


