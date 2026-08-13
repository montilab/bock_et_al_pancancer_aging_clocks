library(SummarizedExperiment)
library(S4Vectors)
library(dplyr)
library(stringr)

# =============================================================================
# CONFIGURATION
# =============================================================================

GISTIC_ROOT <- "/restricted/projectnb/agedisease/projects/pancancer_aging/TCGA-GDC/scna/GISTIC"
TCGA_ROOT <- file.path(
  Sys.getenv("CBMDATA", "/restricted/projectnb/agedisease/CBMrepositoryData/"),
  "TCGA-GDC"
)

path_raw_data <- file.path(GISTIC_ROOT, "raw_data")
path_processed_SE <- file.path(GISTIC_ROOT, "processed_data")
path_metadata <- file.path(TCGA_ROOT, "metadata/data/harmonized_metadata.rds")
path_protein_genes <- file.path(TCGA_ROOT, "metadata/resources/version_115_ensembl_protein_coding_autosomal_genes.RDS")

dir.create(path_processed_SE, recursive = TRUE, showWarnings = FALSE)

projects <- c(
  "TCGA-BRCA", "TCGA-OV", "TCGA-GBM", "TCGA-LUAD", "TCGA-UCEC",
  "TCGA-KIRC", "TCGA-HNSC", "TCGA-LGG", "TCGA-THCA", "TCGA-LUSC",
  "TCGA-PRAD", "TCGA-SKCM", "TCGA-COAD", "TCGA-STAD", "TCGA-BLCA",
  "TCGA-LIHC", "TCGA-CESC", "TCGA-KIRP", "TCGA-SARC", "TCGA-LAML",
  "TCGA-ESCA", "TCGA-PAAD", "TCGA-PCPG", "TCGA-READ", "TCGA-TGCT",
  "TCGA-THYM", "TCGA-KICH", "TCGA-ACC", "TCGA-MESO", "TCGA-UVM",
  "TCGA-DLBC", "TCGA-UCS", "TCGA-CHOL"
)

default_threshold <- 30
cancer_thresholds <- list(
  "TCGA-THYM" = 15,
  "TCGA-THCA" = 8,
  "TCGA-LAML" = 8
)

get_threshold <- function(project) {
  if (project %in% names(cancer_thresholds)) {
    return(cancer_thresholds[[project]])
  }
  return(default_threshold)
}

get_primary_sample_filter <- function(barcodes, project) {
  sample_types <- substr(barcodes, 14, 15)
  if (project == "TCGA-LAML") {
    return(sample_types == "03")
  } else {
    return(sample_types == "01")
  }
}

# =============================================================================
# LOAD REFERENCE DATA
# =============================================================================

cat("==========================================================\n")
cat(" LOADING REFERENCE DATA\n")
cat("==========================================================\n")

if (!file.exists(path_metadata)) {
  stop("Metadata file not found at: ", path_metadata)
}
cat(" -> Loading Harmonized Metadata...\n")
metadata <- readRDS(path_metadata)
metadata$grade <- metadata$histological_grade

if (!file.exists(path_protein_genes)) {
  stop("Protein-coding genes file not found at: ", path_protein_genes)
}
cat(" -> Loading Protein-Coding Autosomal Genes Reference...\n")
autosomal_protein_genes <- readRDS(path_protein_genes)
autosomal_protein_genes <- autosomal_protein_genes %>%
  dplyr::filter(hgnc_symbol != "")

cat(" -> Reference data loaded successfully.\n\n")

# =============================================================================
# PART 1: SEGMENT-LEVEL SE CREATION
# =============================================================================

cat("==========================================================\n")
cat(" PART 1: SEGMENT-LEVEL SE CREATION\n")
cat("==========================================================\n\n")

segment_summary <- data.frame()

for (project in projects) {
  cat(paste0("\n----------------------------------------------------------\n"))
  cat(paste0(" PROCESSING SEGMENT-LEVEL: ", project, "\n"))

  project_dir <- file.path(path_raw_data, project)
  all_lesions_path <- file.path(project_dir, "all_lesions.conf_99.txt")
  amp_genes_path <- file.path(project_dir, "table_amp.conf_99.txt")
  del_genes_path <- file.path(project_dir, "table_del.conf_99.txt")

  if (!file.exists(all_lesions_path)) {
    cat(" -> GISTIC files not found - skipping.\n")
    next
  }

  tryCatch({
    mat_raw <- read.delim(all_lesions_path, stringsAsFactors = FALSE)
    table_amp <- read.delim(amp_genes_path, stringsAsFactors = FALSE)
    table_del <- read.delim(del_genes_path, stringsAsFactors = FALSE)

    n_amps <- nrow(table_amp)
    n_dels <- nrow(table_del)
    total_regions <- n_amps + n_dels

    mat_discrete_raw <- mat_raw[seq_len(total_regions), ]
    mat_cont_raw <- mat_raw[(total_regions + 1):nrow(mat_raw), ]

    if (n_amps > 0) {
      mat_amp_disc <- mat_discrete_raw[seq_len(n_amps), ]
      mat_amp_disc$Unique.Name <- gsub("\\s+", "_", mat_amp_disc$Unique.Name)
      mat_amp_disc$index <- paste0(mat_amp_disc$Unique.Name, "_", mat_amp_disc$Descriptor)
      mat_amp_combined <- cbind(mat_amp_disc %>% dplyr::select(index), table_amp)
    }

    if (n_dels > 0) {
      mat_del_disc <- mat_discrete_raw[(n_amps + 1):total_regions, ]
      mat_del_disc$Unique.Name <- gsub("\\s+", "_", mat_del_disc$Unique.Name)
      mat_del_disc$index <- paste0(mat_del_disc$Unique.Name, "_", mat_del_disc$Descriptor)
      mat_del_combined <- cbind(mat_del_disc %>% dplyr::select(index), table_del)
    }

    if (n_amps > 0 && n_dels > 0) {
      rowdata_full <- rbind(mat_amp_combined, mat_del_combined)
    } else if (n_amps > 0 && n_dels == 0) {
      rowdata_full <- mat_amp_combined
    } else if (n_dels > 0 && n_amps == 0) {
      rowdata_full <- mat_del_combined
    } else {
      stop("GISTIC detected no significant peaks for this cancer type")
    }

    rowdata_full$index <- gsub("\\s+", "", rowdata_full$index)
    rownames(rowdata_full) <- rowdata_full$index

    rowdata_full$genes_in_region <- rowdata_full$genes_in_region %>%
      str_replace_all("\\|ENSG[^,]*|,$", "") %>%
      str_remove_all("[\\[\\]]")
    rowdata_full$genes_in_peak <- rowdata_full$genes_in_peak %>%
      str_replace_all("\\|ENSG[^,]*|,$", "") %>%
      str_remove_all("[\\[\\]]")

    process_matrix_cols <- function(mat_subset, proj) {
      data_cols <- mat_subset[, 10:ncol(mat_subset), drop = FALSE]
      colnames(data_cols) <- gsub("\\.", "-", colnames(data_cols))
      keep <- get_primary_sample_filter(colnames(data_cols), proj)
      data_cols <- data_cols[, keep, drop = FALSE]
      colnames(data_cols) <- substr(colnames(data_cols), 1, 12)
      if (any(duplicated(colnames(data_cols)))) {
        data_cols <- data_cols[, !duplicated(colnames(data_cols)), drop = FALSE]
      }
      as.matrix(data_cols)
    }

    mat_3_level <- process_matrix_cols(mat_discrete_raw, project)
    rownames(mat_3_level) <- rowdata_full$index

    mat_cont <- process_matrix_cols(mat_cont_raw, project)
    rownames(mat_cont) <- rowdata_full$index

    mat_2_level <- mat_3_level
    mat_2_level[mat_2_level == 2] <- 1

    common_local <- intersect(colnames(mat_3_level), colnames(mat_cont))
    mat_3_level <- mat_3_level[, common_local, drop = FALSE]
    mat_2_level <- mat_2_level[, common_local, drop = FALSE]
    mat_cont <- mat_cont[, common_local, drop = FALSE]

    common_meta <- intersect(colnames(mat_3_level), rownames(metadata))

    if (length(common_meta) == 0) {
      cat(" -> No samples match harmonized metadata - skipping.\n")
      next
    }

    project_specific_metadata <- metadata[common_meta, , drop = FALSE]
    project_specific_metadata <- DataFrame(project_specific_metadata)

    if ("age_at_diagnosis" %in% colnames(project_specific_metadata)) {
      project_specific_metadata$age_at_diagnosis_years <- project_specific_metadata$age_at_diagnosis / 365
    }

    target_order <- rownames(project_specific_metadata)
    mat_3_level <- mat_3_level[, target_order, drop = FALSE]
    mat_2_level <- mat_2_level[, target_order, drop = FALSE]
    mat_cont <- mat_cont[, target_order, drop = FALSE]

    se <- SummarizedExperiment::SummarizedExperiment(
      assays = list(
        GISTIC_thresholded_three_level = mat_3_level,
        GISTIC_thresholded_two_level = mat_2_level,
        GISTIC_continuous = mat_cont
      ),
      rowData = rowdata_full,
      colData = project_specific_metadata
    )

    output_filename <- paste0(project, "_SCNA_GISTIC_segment.rds")
    saveRDS(se, file = file.path(path_processed_SE, output_filename))

    segment_summary <- rbind(segment_summary, data.frame(
      project = project,
      n_amp_peaks = n_amps,
      n_del_peaks = n_dels,
      n_samples = ncol(se),
      status = "Success"
    ))

    cat(paste0(
      " -> Saved: ", output_filename, " | Samples: ", ncol(se),
      " | Amp Peaks: ", n_amps, " | Del Peaks: ", n_dels, "\n"
    ))
  }, error = function(e) {
    segment_summary <<- rbind(segment_summary, data.frame(
      project = project,
      n_amp_peaks = NA,
      n_del_peaks = NA,
      n_samples = NA,
      status = paste0("ERROR: ", e$message)
    ))
    cat(paste0(" -> ERROR: ", e$message, "\n"))
  })
}

cat("\n Segment-level processing complete.\n\n")

# =============================================================================
# PART 2: GENE-LEVEL SE CREATION
# =============================================================================

cat("==========================================================\n")
cat(" PART 2: GENE-LEVEL SE CREATION\n")
cat("==========================================================\n\n")

gene_summary <- data.frame()

for (project in projects) {
  cat(paste0("\n----------------------------------------------------------\n"))
  cat(paste0(" PROCESSING GENE-LEVEL: ", project, "\n"))

  input_file <- file.path(path_raw_data, project, "all_thresholded.by_genes.txt")

  if (!file.exists(input_file)) {
    cat(" -> Gene-level file not found - skipping.\n")
    next
  }

  tryCatch({
    mat_raw <- read.delim(input_file, stringsAsFactors = FALSE)

    rowdata_df <- mat_raw[, 1:3]
    rownames(rowdata_df) <- rowdata_df$Gene.Symbol

    mat <- as.matrix(mat_raw[, !colnames(mat_raw) %in% c("Gene.Symbol", "Locus.ID", "Cytoband")])
    rownames(mat) <- rowdata_df$Gene.Symbol
    colnames(mat) <- gsub("\\.", "-", colnames(mat))

    keep_cols <- get_primary_sample_filter(colnames(mat), project)
    if (sum(keep_cols) == 0) {
      cat(" -> No primary tumor samples found.\n")
      next
    }

    mat <- mat[, keep_cols, drop = FALSE]
    colnames(mat) <- substr(colnames(mat), 1, 12)

    if (any(duplicated(colnames(mat)))) {
      mat <- mat[, !duplicated(colnames(mat)), drop = FALSE]
    }

    common_samples <- intersect(colnames(mat), rownames(metadata))

    if (length(common_samples) == 0) {
      cat(" -> No samples match harmonized metadata - skipping.\n")
      next
    }

    mat <- mat[, common_samples, drop = FALSE]
    project_specific_metadata <- metadata[common_samples, , drop = FALSE]
    project_specific_metadata <- DataFrame(project_specific_metadata)

    if ("age_at_diagnosis" %in% colnames(project_specific_metadata)) {
      project_specific_metadata$age_at_diagnosis_years <- project_specific_metadata$age_at_diagnosis / 365
    }

    mat <- mat[, match(rownames(project_specific_metadata), colnames(mat)), drop = FALSE]

    genes_before <- nrow(mat)
    mat <- mat[rownames(mat) %in% autosomal_protein_genes$hgnc_symbol, , drop = FALSE]
    rowdata_df <- rowdata_df[rownames(mat), , drop = FALSE]
    genes_after <- nrow(mat)

    if (genes_after == 0) {
      cat(" -> No genes remain after protein-coding filter - skipping.\n")
      next
    }

    freq_threshold <- get_threshold(project)
    cat(paste0(" -> Using frequency threshold: ", freq_threshold, "%\n"))

    total_samples <- ncol(mat)
    cna_levels <- c("-2", "-1", "0", "1", "2")
    cna_counts <- t(apply(mat, 1, function(x) table(factor(x, levels = cna_levels))))

    loss_pct <- (cna_counts[, "-2"] + cna_counts[, "-1"]) / total_samples * 100
    gain_pct <- (cna_counts[, "2"] + cna_counts[, "1"]) / total_samples * 100

    amplified_genes <- rownames(cna_counts)[loss_pct <= 10 & gain_pct >= freq_threshold]
    deleted_genes <- rownames(cna_counts)[gain_pct <= 10 & loss_pct >= freq_threshold]

    genes_to_keep <- c(amplified_genes, deleted_genes)

    if (length(genes_to_keep) == 0) {
      cat(paste0(
        " -> No genes pass frequency filtering (Amp: ", length(amplified_genes),
        ", Del: ", length(deleted_genes), ") - skipping.\n"
      ))

      gene_summary <- rbind(gene_summary, data.frame(
        project = project,
        threshold_used = freq_threshold,
        genes_before_filter = genes_before,
        genes_after_pc_filter = genes_after,
        n_amplified = 0,
        n_deleted = 0,
        n_samples = ncol(mat),
        status = "No genes pass filter"
      ))
      next
    }

    mat <- mat[genes_to_keep, , drop = FALSE]
    rowdata_df <- rowdata_df[rownames(mat), , drop = FALSE]

    mat_five <- mat

    mat_amp_3 <- mat_five
    mat_amp_3[mat_amp_3 <= 0] <- 0
    mat_del_3 <- mat_five
    mat_del_3[mat_del_3 >= 0] <- 0
    mat_del_3 <- abs(mat_del_3)

    mat_amp_2 <- mat_amp_3
    mat_amp_2[mat_amp_2 == 2] <- 1
    mat_del_2 <- mat_del_3
    mat_del_2[mat_del_2 == 2] <- 1

    se <- SummarizedExperiment::SummarizedExperiment(
      assays = list(
        five_levels = mat_five,
        amplified_three_levels = mat_amp_3,
        deleted_three_levels = mat_del_3,
        amplified_two_levels = mat_amp_2,
        deleted_two_levels = mat_del_2
      ),
      rowData = rowdata_df,
      colData = project_specific_metadata
    )

    metadata(se)$frequency_threshold <- freq_threshold

    output_filename <- paste0(project, "_SCNA_GISTIC_gene.rds")
    saveRDS(se, file = file.path(path_processed_SE, output_filename))

    gene_summary <- rbind(gene_summary, data.frame(
      project = project,
      threshold_used = freq_threshold,
      genes_before_filter = genes_before,
      genes_after_pc_filter = genes_after,
      n_amplified = length(amplified_genes),
      n_deleted = length(deleted_genes),
      n_samples = ncol(se),
      status = "Success"
    ))

    cat(paste0(
      " -> Saved: ", output_filename, " | Genes: ", nrow(se),
      " (Amp: ", length(amplified_genes), ", Del: ", length(deleted_genes),
      ") | Samples: ", ncol(se), " | Threshold: ", freq_threshold, "%\n"
    ))
  }, error = function(e) {
    gene_summary <<- rbind(gene_summary, data.frame(
      project = project,
      threshold_used = get_threshold(project),
      genes_before_filter = NA,
      genes_after_pc_filter = NA,
      n_amplified = NA,
      n_deleted = NA,
      n_samples = NA,
      status = paste0("ERROR: ", e$message)
    ))
    cat(paste0(" -> ERROR: ", e$message, "\n"))
  })
}

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n==========================================================\n")
cat(" PROCESSING COMPLETE\n")
cat("==========================================================\n\n")

cat("SEGMENT-LEVEL SUMMARY:\n")
print(segment_summary)

cat("\nGENE-LEVEL SUMMARY:\n")
print(gene_summary)

cat("\n\nOutput files saved to:", path_processed_SE, "\n")
cat("\nNEXT STEP:\n")
cat("Payton runs 03_Covariates.Rmd to generate SCNA_covariates.rds\n")
cat("Update se_directory to: file.path(TCGA, \"scna/GISTIC/processed_data/\")\n")

cat("\n==========================================================\n")
cat(" ALL JOBS COMPLETED\n")
cat("==========================================================\n")