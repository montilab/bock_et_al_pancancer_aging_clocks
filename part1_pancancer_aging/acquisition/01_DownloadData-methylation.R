# Author: Payton Bock
## Downloading TCGA methylation data using TCGAbiolinks
## This version of the script is built for qsub job submission using 01_DownloadData_JobSubmission.sh
# Run once per project via 01_DownloadData_JobSubmission.sh, across all 33 TCGA projects:
# TCGA-BRCA, TCGA-OV, TCGA-GBM, ... (same list as the RNA-seq and miRNA scripts)

# Get arguments from the command line
args <- commandArgs(trailingOnly = TRUE)

# Check if the argument is provided
if (length(args) == 0) {
  stop("No project name provided. Please provide a project name as an argument.")
}

# Use the argument as the project name
project_name <- args[1]

# Load TCGA biolinks library
.libPaths("/restricted/projectnb/agedisease/projects/pancancer_aging_pbock/r-library/")
library(TCGAbiolinks)
library(sesameData)
library(sesame)

# Determines paths for outputs
# Sets the CBMDATA directory path
Sys.setenv(CBMDATA = "/restricted/projectnb/agedisease/projects/pancancer_aging")

# Combines CBMDATA path with sub-directories for data output
PATH <- file.path(Sys.getenv("CBMDATA"), "TCGA-GDC", "methylation", "raw_data")
PATH_rawdl <- file.path(PATH, "project_downloads")
PATH_se <- file.path(PATH, "raw_SE")

# Function to download data from TCGA
# Takes in project name from input to query
# Platform and Data type is kept consistent. Adjust if necessary
download_tcga_data <- function(project) {
  if (project == "TCGA-LAML") {
    query <- TCGAbiolinks::GDCquery(
      project = "TCGA-LAML",
      data.category = "DNA Methylation",
      platform = "Illumina Human Methylation 450",
      data.type = "Methylation Beta Value"
    )
  } else {
    query <- TCGAbiolinks::GDCquery(
      project = project,
      data.category = "DNA Methylation",
      platform = "Illumina Human Methylation 450",
      data.type = "Methylation Beta Value",
      sample.type = "Primary Tumor"
    )
  }
  message("Starting download for project: ", project_name)
  
  # Error handling for download
  tryCatch({
    TCGAbiolinks::GDCdownload(query, directory = PATH_rawdl)
    data <- TCGAbiolinks::GDCprepare(query, directory = PATH_rawdl)
  }, error = function(e) {
    message("An error occurred: ", e$message)
    stop("Failed to download or prepare data.")
  })
  
  return(data)
}

# Creates SE object using data (after first downloading the data)
tcga_data <- download_tcga_data(project_name)

# Prepares name of summarizedexperiment rds file
SE_output_filename <- paste0(project_name, "_methylation.rds")

# Save SummarizedExperiment to proper folder
saveRDS(tcga_data, file = file.path(PATH_se, SE_output_filename))
