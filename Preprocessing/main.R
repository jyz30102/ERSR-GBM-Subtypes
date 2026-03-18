#!/usr/bin/env Rscript
#' =============================================================================
#' Project: GBM/LGG Methylation Preprocessing Pipeline
#' Module:  main.R
#' Purpose: Orchestrate the complete preprocessing workflow for GSE90496
#' =============================================================================

# :: Load dependencies ---------------------------------------------------------
suppressPackageStartupMessages({
  library(GEOquery)
  library(minfi)
  library(limma)
  library(yaml)
})

# :: Source project modules ----------------------------------------------------
source("R/config.R")
source("R/utils.R")

# :: Load configuration --------------------------------------------------------
config <- .load_config("config/config.yaml")
config <- .init_logging(config)

.log_msg("=== Starting Methylation Preprocessing Pipeline ===", "INFO")
.log_msg(sprintf("Project: %s v%s", config$project$name, config$project$version))

# :: Step 1: Download and parse sample annotation ------------------------------
.log_msg("Step 1: Downloading sample annotation from GEO...", "INFO")

gse <- getGEO(
  config$geo$accession, 
  GSEMatrix = TRUE, 
  getGPL = FALSE
)

anno <- pData(gse[[1]])  # Extract phenotype data
.log_msg(sprintf("Downloaded annotation: %d samples", nrow(anno)), "INFO")

# Filter for target classes (LGG + GBM)
target_classes <- config$geo$target_classes
class_col <- config$geo$methylation_class_col

if (!class_col %in% colnames(anno)) {
  stop(sprintf("Column '%s' not found in annotation. Available: %s", 
               class_col, paste(colnames(anno), collapse = ", ")))
}

anno_filtered <- anno[anno[[class_col]] %in% target_classes, , drop = FALSE]
.log_msg(sprintf("Filtered to target classes %s: %d samples", 
                 paste(target_classes, collapse = "/"), nrow(anno_filtered)), "INFO")

# :: Step 2: Read raw IDAT files -----------------------------------------------
.log_msg("Step 2: Reading raw IDAT files...", "INFO")

# Parse file paths from GEO supplementary files
# Format: extract filename from supplementary_file column
extract_idat_path <- function(supp_file, raw_dir) {
  # Example: "https://.../GSE90496_RAW/.../sample_Grn.idat"
  # We want: "data/raw/GSE90496_RAW/sample"
  fname <- basename(supp_file)
  sample_name <- gsub("_Grn.*", "", fname)
  file.path(raw_dir, sample_name)
}

idat_paths <- sapply(
  anno_filtered$supplementary_file, 
  extract_idat_path, 
  raw_dir = config$paths$raw_idat,
  USE.NAMES = FALSE
)

# Verify files exist
missing_files <- idat_paths[!file.exists(paste0(idat_paths, "_Grn.idat"))]
if (length(missing_files) > 0) {
  stop(sprintf("Missing IDAT files: %d samples. First: %s", 
               length(missing_files), missing_files[1]))
}

# Read IDAT files using minfi
start_time <- Sys.time()
RGset <- read.metharray(idat_paths, verbose = TRUE)
elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
.log_msg(sprintf("Read %d samples in %.1f minutes", ncol(RGset), elapsed), "INFO")

# :: Step 3: Illumina normalization --------------------------------------------
.log_msg("Step 3: Running Illumina normalization...", "INFO")

# Source custom preprocessing functions
if (file.exists(config$paths$functions_script)) {
  source(config$paths$functions_script)
  .log_msg(sprintf("Loaded custom functions: %s", config$paths$functions_script))
} else {
  warning(sprintf("Custom functions script not found: %s", config$paths$functions_script))
}

# Run normalization (assumes MNPpreprocessIllumina is defined in functions script)
start_time <- Sys.time()
Mset <- MNPpreprocessIllumina(RGset)
elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
.log_msg(sprintf("Normalization complete in %.1f minutes", elapsed), "INFO")

# :: Step 4: Probe filtering ---------------------------------------------------
.log_msg("Step 4: Filtering probes...", "INFO")

# Load filter lists
filters <- .load_filter_probes(
  config$paths$filter_dir, 
  config$preprocessing$filters
)

# Get indices to remove
remove_idx <- .get_remove_indices(
  rownames(Mset), 
  filters, 
  config$preprocessing$filters
)

.log_msg(sprintf("Total probes to remove: %d / %d", 
                 length(remove_idx), nrow(Mset)), "INFO")

# Apply filtering
Mset_filtered <- Mset[-remove_idx, ]
.log_msg(sprintf("Retained probes: %d", nrow(Mset_filtered)), "INFO")

# Save intermediate result
if (config$output$save_intermediate) {
  .save_rdata(
    list(Mset_filtered = Mset_filtered, anno = anno_filtered),
    file.path(config$paths$output_dir, "Mset_filtered.RData"),
    config$output$compression
  )
}

# Clean up memory
rm(Mset, RGset); gc()

# :: Step 5: Batch effect correction -------------------------------------------
.log_msg("Step 5: Correcting batch effects...", "INFO")

# Extract methylation/unmethylation signals
methy <- getMeth(Mset_filtered)
unmethy <- getUnmeth(Mset_filtered)
rm(Mset_filtered); gc()

# Define batch variable from annotation
batch_col <- config$preprocessing$batch_correction$batch_col
batch_levels <- config$preprocessing$batch_correction$batch_levels

if (!batch_col %in% colnames(anno_filtered)) {
  stop(sprintf("Batch column '%s' not found in annotation", batch_col))
}

# Create numeric batch vector
batch_vec <- anno_filtered[[batch_col]]
batch_numeric <- recode(batch_vec, !!!batch_levels)
batch_numeric <- as.numeric(as.factor(batch_numeric))  # Ensure numeric

.log_msg(sprintf("Batch distribution: %s", 
                 paste(table(batch_vec), collapse = ", ")), "INFO")

# Apply limma's removeBatchEffect
bc_cfg <- config$preprocessing$batch_correction
pseudocount <- bc_cfg$pseudocount

start_time <- Sys.time()
methy_log <- log2(methy + pseudocount)
unmethy_log <- log2(unmethy + pseudocount)

methy_bc <- 2^removeBatchEffect(methy_log, batch = batch_numeric)
unmethy_bc <- 2^removeBatchEffect(unmethy_log, batch = batch_numeric)

elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
.log_msg(sprintf("Batch correction complete in %.1f minutes", elapsed), "INFO")

# Extract batch effect coefficients (for adjusting new samples)
# Use first sample of each batch as reference
batch_unique <- unique(batch_vec)
coef_list <- list(methy = list(), unmethy = list())

for (b in batch_unique) {
  idx <- which(batch_vec == b)[1]  # First sample of this batch
  coef_list$methy[[b]] <- log2(methy_bc[, idx]) - log2(methy[, idx] + pseudocount)
  coef_list$unmethy[[b]] <- log2(unmethy_bc[, idx]) - log2(unmethy[, idx] + pseudocount)
}

# Save batch coefficients
if (config$output$save_intermediate) {
  .save_rdata(
    coef_list,
    file.path(config$paths$output_dir, "ba.coef.RData"),
    config$output$compression
  )
}

# :: Step 6: Calculate beta values ---------------------------------------------
.log_msg("Step 6: Calculating beta values...", "INFO")

beta_cfg <- config$preprocessing$beta
betas <- .calc_beta_values(methy_bc, unmethy_bc, beta_cfg$offset)

# Transpose to samples × probes (standard for downstream analysis)
betas_df <- as.data.frame(t(betas))
rownames(betas_df) <- rownames(anno_filtered)  # Ensure sample IDs match

.log_msg(sprintf("Beta matrix: %d samples × %d probes", 
                 nrow(betas_df), ncol(betas_df)), "INFO")

# :: Step 7: Save final results ------------------------------------------------
.log_msg("Step 7: Saving final results...", "INFO")

output_dir <- config$paths$output_dir

# Save as RData
if (config$output$save_final_beta) {
  .save_rdata(
    list(betas = betas_df, anno = anno_filtered),
    file.path(output_dir, "betas.ba.RData"),
    config$output$compression
  )
  
  # Optional: Save as CSV for interoperability
  if (config$output$file_format == "CSV" || TRUE) {  # Always save CSV for flexibility
    .save_beta_csv(
      betas_df, 
      file.path(output_dir, "betas.ba.csv"),
      anno_filtered
    )
  }
}

# Save processing summary
summary <- list(
  timestamp = Sys.time(),
  config_version = config$project$version,
  input_samples = nrow(anno),
  filtered_samples = nrow(anno_filtered),
  input_probes = nrow(Mset),
  filtered_probes = nrow(Mset_filtered),
  probes_removed = length(remove_idx),
  batch_levels = batch_levels,
  beta_offset = beta_cfg$offset,
  output_files = list(
    rdata = file.path(output_dir, "betas.ba.RData"),
    csv = file.path(output_dir, "betas.ba.csv"),
    coef = file.path(output_dir, "ba.coef.RData")
  )
)

saveRDS(summary, file.path(output_dir, "processing_summary.rds"))
.log_msg("Pipeline completed successfully!", "INFO")

# :: Return results (for interactive use) --------------------------------------
invisible(list(
  betas = betas_df,
  annotation = anno_filtered,
  batch_coef = coef_list,
  summary = summary
))
