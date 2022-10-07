library(here)

runs <- c(
  "billi2022_seurat-object", "BF_LE_08_GD_pre_h"
)

# rmarkdown::render(
#   "analysis-report-samples.Rmd",
#   output_file = here("figs", "my_markdown"),
#   params = list(run = runs[1])
# )

filenames <- paste0(runs, "-report.html")
params <- map(runs, ~list(run = .))

for(i in seq_along(filenames)){
  
  quarto::quarto_render(here("reports","2022-09-05_analysis-report-samples.qmd"), output_file = here("reports", filenames[i]), execute_params = params[[i]])
  
}

map2(
  filenames,
  params,
  ~quarto::quarto_render(here("reports","2022-09-05_analysis-report-samples.qmd"), output_file = here("reports", .x), execute_params = .y )
)

# quarto::quarto_render(here("reports","2022-09-05_analysis-report-samples.qmd"), output_file = here("figs", filenames[2]), execute_params = params[[2]], as_job = "quarto.render_as_job" )

