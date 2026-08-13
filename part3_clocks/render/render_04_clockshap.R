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
  
  message("Rendering 04_clockshap: ", tcga, " / ", gtex)
  
  rmarkdown::render(
    input = "part3_clocks/Rmd/04_clockshap_template.Rmd",
    params = list(
      tcga_project = tcga,
      gtex_project = gtex
    ),
    output_file = file.path(
      "..", "html",
      paste0("04_clockshap_", tcga, "_", gtex, ".html")
    ),
    clean = TRUE,
    envir = new.env(),
    output_options = list(
      css = here::here("R", "style.css")
    )
  )
}
