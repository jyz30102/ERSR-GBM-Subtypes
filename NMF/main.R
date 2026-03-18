#!/usr/bin/env Rscript
#' =============================================================================
#' Project: ERSR-GBM-Subtypes
#' Module:  NMF/main.R
#' Purpose: NMF Analysis Entry Point
#' =============================================================================

suppressPackageStartupMessages({
  library(NMF)
  library(tinyarray)
  library(ggplot2)
  library(dplyr)
  library(Cairo)
  library(yaml)
})

# :: Load modules --------------------------------------------------------------
source("NMF/utils.R")
source("NMF/nmf_analysis.R")

# :: Load configuration --------------------------------------------------------
config <- .load_config("NMF/config.yaml")
.log_msg("=== Starting NMF Analysis Pipeline ===", "INFO")

# :: Step 1: Load preprocessed data --------------------------------------------
.log_msg("Step 1: Loading preprocessed beta matrix...", "INFO")
beta_matrix <- .load_beta_matrix(config$paths$input_data)
beta_matrix <- .validate_beta_matrix(beta_matrix)

# :: Step 2: Rank selection ----------------------------------------------------
.log_msg("Step 2: Performing rank selection...", "INFO")

rank_range <- seq(
  config$nmf$rank_range$min,
  config$nmf$rank_range$max,
  by = config$nmf$rank_range$step
)

rank_results <- .nmf_rank_selection(
  beta_matrix = beta_matrix,
  rank_range = rank_range,
  nrun = config$nmf$nrun,
  seed = config$nmf$seed,
  method = config$nmf$method
)

# :: Step 3: Visualize rank selection ------------------------------------------
.log_msg("Step 3: Generating rank selection plots...", "INFO")

.plot_cophenetic(
  rank_results,
  file.path(config$paths$figures_dir, "cophenetic_rank_selection"),
  config
)

# :: Step 4: Final NMF run -----------------------------------------------------
.log_msg("Step 4: Running final NMF with optimal rank...", "INFO")

optimal_rank <- rank_results$optimal_rank
nmf_fit <- .nmf_final_run(
  beta_matrix = beta_matrix,
  rank = optimal_rank,
  nrun = config$nmf$nrun,
  seed = config$nmf$seed,
  method = config$nmf$method
)

# :: Step 5: Get sample assignments --------------------------------------------
.log_msg("Step 5: Extracting sample assignments...", "INFO")
assignments <- .get_sample_assignments(nmf_fit)

# :: Step 6: Visualize results -------------------------------------------------
.log_msg("Step 6: Generating visualization...", "INFO")

.plot_consensus(
  nmf_fit,
  file.path(config$paths$figures_dir, "consensus_matrix"),
  config
)

.plot_clusters(
  assignments,
  anno = NULL,
  file.path(config$paths$figures_dir, "cluster_distribution"),
  config
)

# :: Step 7: Save results ------------------------------------------------------
.log_msg("Step 7: Saving results...", "INFO")

.save_nmf_results(
  nmf_fit = nmf_fit,
  assignments = assignments,
  output_dir = config$paths$output_dir,
  config = config
)

# :: Complete ------------------------------------------------------------------
.log_msg("=== NMF Analysis Completed Successfully ===", "INFO")
.log_msg(sprintf("Optimal rank: k = %d", optimal_rank))
.log_msg(sprintf("Cluster distribution: %s", 
                 paste(table(assignments), collapse = ", ")))
.log_msg(sprintf("Results saved to: %s", config$paths$output_dir))
