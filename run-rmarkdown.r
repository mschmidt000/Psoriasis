library(here)

runs <- c(
  "BF-LE-01-KT-PSO_all",  "BF-LE-02-PG-PSO_all",  "BF_LE_03_VC_03_all", "BF_LE_06_KS_LE_all",  "BF_LE_08_GD_pre_h"
)

# rmarkdown::render(
#   "analysis-report-samples.Rmd",
#   output_file = here("figs", "my_markdown"),
#   params = list(run = runs[1])
# )

filenames <- paste0(runs, "-report.html")
params <- map(runs, ~list(run = .))

map2(
  filenames,
  params,
  ~rmarkdown::render("analysis-report-samples.Rmd", output_file = here::here("figs", .x), params = .y )
)
