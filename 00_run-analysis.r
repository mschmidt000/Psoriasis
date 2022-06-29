### run analysis of melanoma scrna-seq data
### 19.05.22
library(here)

source(here("src", "01_load-data-and-preprocess.r"))
source(here("src", "02_reduce-dimension-and-cluster.r"))
source(here("src", "03_predict-cell-type.r"))

source(here("src", "04_integrate-data.r"))
source(here("src", "05_subclustering-melanocytes.r"))
source(here("src", "06_receptor-ligand-analysis.r"))
source(here("src", "06_run-opossom.r"))
