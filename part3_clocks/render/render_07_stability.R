map <- list(
  BRCA = "BREAST",
  LUAD = "LUNG",
  LUSC = "LUNG",
  COAD = "COLON",
  LGG  = "BRAIN",
  GBM  = "BRAIN",
  KIRC = "KIDNEY",
  HNSC = "ESOPHAGUS"
)

jobs <- data.frame(
  tcga_project = names(map),
  gtex_project = unlist(map),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(jobs))) {
  tcga <- jobs$tcga_project[i]
  gtex <- jobs$gtex_project[i]
  
  message("Rendering stability: ", tcga, " / ", gtex)
  
  rmarkdown::render(
    input = "part3_clocks/Rmd/07_stability_clean_template.Rmd",
    params = list(
      tcga_project = tcga,
      gtex_project = gtex,
      feature_set  = "PC",
      out_root     = "results/derived/GTEx",
      min_nat      = 10
    ),
    output_file = file.path(
      "..", "html",
      paste0("07_stability_", tcga, "_", gtex, ".html")
    ),
    clean = TRUE,
    envir = new.env(),
    output_options = list(
      css = here::here("R", "style.css")
    )
  )
}