# This script renders the Mutation EDA & Filtering reports for all cancer types
# Location: /restricted/projectnb/agedisease/projects/pancancer_aging/scripts/mutation/

library(tidyverse)
library(rmarkdown)

DATA_ROOT <- "/restricted/projectnb/agedisease/projects/pancancer_aging/TCGA-GDC/mutation"
SCRIPT_ROOT <- "/restricted/projectnb/agedisease/projects/pancancer_aging/scripts/mutation"
RESULTS_ROOT <- "/restricted/projectnb/agedisease/projects/pancancer_aging/results/mutation"

render_reports <- function(rmd, output_dir, covariate_df, suffix) {
  required_params <- c(
    "cancer_type", "include_gender", "include_race", "races_included",
    "include_subtype", "include_HPV_status", "include_grade",
    "grades_included", "include_purity"
  )
  missing <- setdiff(required_params, names(covariate_df))
  if (length(missing) > 0) stop("Missing columns in covariate_df: ", paste(missing, collapse = ", "))

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  for (i in seq_len(nrow(covariate_df))) {
    params_list <- as.list(covariate_df[i, ])
    output_file <- file.path(output_dir, paste0(params_list$cancer_type, "_", suffix))

    print(paste0("[", i, "/", nrow(covariate_df), "] Rendering: ", params_list$cancer_type, " -> ", output_file))

    tryCatch(
      {
        rmarkdown::render(
          input = rmd,
          output_file = output_file,
          params = params_list,
          quiet = TRUE
        )
        print("Done!")
      },
      error = function(e) {
        message(paste0("ERROR rendering ", params_list$cancer_type, ": ", e$message))
      }
    )
  }
}

covariates_mutation <- readRDS(file.path(
  "/restricted/projectnb/agedisease/projects/pancancer_aging",
  "TCGA-GDC", "metadata", "covariates", "MUTATION_covariates.rds"
)) %>%
  dplyr::select(-omics_type, -num_samples) %>%
  mutate(
    include_race = if_else(cancer_type %in% c("TCGA-COAD", "TCGA-READ"), FALSE, include_race)
  )

render_reports(
  rmd = file.path(SCRIPT_ROOT, "03_EDA_and_filtering-mutation.Rmd"),
  output_dir = file.path(RESULTS_ROOT, "htmls"),
  covariate_df = covariates_mutation,
  suffix = "Mutation_EDA_and_filtering.html"
)
