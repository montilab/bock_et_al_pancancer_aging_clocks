library(TCGAbiolinks)
library(stringr)

projects <- c(
  "TCGA-BRCA", "TCGA-OV", "TCGA-GBM", "TCGA-LUAD", "TCGA-UCEC",
  "TCGA-KIRC", "TCGA-HNSC", "TCGA-LGG", "TCGA-THCA", "TCGA-LUSC",
  "TCGA-PRAD", "TCGA-SKCM", "TCGA-COAD", "TCGA-STAD", "TCGA-BLCA",
  "TCGA-LIHC", "TCGA-CESC", "TCGA-KIRP", "TCGA-SARC", "TCGA-LAML",
  "TCGA-ESCA", "TCGA-PAAD", "TCGA-PCPG", "TCGA-READ", "TCGA-TGCT",
  "TCGA-THYM", "TCGA-KICH", "TCGA-ACC", "TCGA-MESO", "TCGA-UVM",
  "TCGA-DLBC", "TCGA-UCS", "TCGA-CHOL"
)

path_scna_root <- "/restricted/projectnb/agedisease/projects/pancancer_aging/TCGA-GDC/scna/GISTIC/raw_data"
dir.create(path_scna_root, recursive = TRUE, showWarnings = FALSE)

original_wd <- getwd()

cat(paste0("Starting GISTIC Download for ", length(projects), " projects.\n"))

for (project_full in projects) {
  project_short <- stringr::str_remove(project_full, "TCGA-")

  cat(paste0("\n==========================================================\n"))
  cat(paste0(" PROCESSING PROJECT: ", project_full, "\n"))

  project_dir <- file.path(path_scna_root, project_full)
  dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)

  cat(" -> Changing directory to:", project_dir, "\n")
  setwd(project_dir)

  cat(" -> Querying Firehose (getGistic)...\n")
  TCGAbiolinks::getGistic(project_short, type = "threshold")

  tar_files <- list.files(pattern = "\\.tar\\.gz$", full.names = TRUE)

  if (length(tar_files) > 0) {
    cat(" -> Found tarball(s). Extracting and flattening...\n")

    for (tgz in tar_files) {
      untar(tgz, exdir = project_dir)

      nested_dirs <- list.files(project_dir,
        pattern = "^gdac\\.broadinstitute",
        full.names = TRUE, include.dirs = TRUE
      )
      nested_dirs <- nested_dirs[dir.exists(nested_dirs)]

      if (length(nested_dirs) > 0) {
        for (nested in nested_dirs) {
          files_inside <- list.files(nested, full.names = TRUE)
          file.rename(from = files_inside, to = file.path(project_dir, basename(files_inside)))
          unlink(nested, recursive = TRUE)
        }
        cat(" -> Folder structure flattened.\n")
      }

      file.remove(tgz)
    }
  }
}

setwd(original_wd)

cat("\n==========================================================\n")
cat(" ALL DOWNLOADS COMPLETED.\n")
cat("==========================================================\n")