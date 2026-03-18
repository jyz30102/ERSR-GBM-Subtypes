#!/usr/bin/env Rscript
#' =============================================================================
#' Module: nmf_analysis.R
#' Purpose: Core NMF analysis functions for GBM subtype discovery
#'          (Includes visualization functions)
#' =============================================================================

suppressPackageStartupMessages({
  library(NMF)
  library(tinyarray)
  library(ggplot2)
  library(dplyr)
  library(Cairo)
  library(yaml)
})

# :: NMF Rank Selection -------------------------------------------------------

.nmf_rank_selection <- function(beta_matrix, rank_range, nrun = 50, 
                                 seed = 42, method = "brunet") {
  .log_msg("Starting NMF rank selection...", "INFO")
  
  cache_file <- file.path("cache/NMF", "rank_selection.RData")
  if (file.exists(cache_file)) {
    .log_msg("Loading cached rank selection results...", "INFO")
    load(cache_file)
    return(rank_results)
  }
  
  results <- list()
  cophenetic_vals <- numeric(length(rank_range))
  dispersion_vals <- numeric(length(rank_range))
  
  for (i in seq_along(rank_range)) {
    k <- rank_range[i]
    .log_msg(sprintf("Testing rank k = %d/%d...", k, max(rank_range)), "INFO")
    
    start_time <- Sys.time()
    
    nmf_fit <- nmf(
      t(beta_matrix),
      rank = k,
      nrun = nrun,
      seed = seed,
      method = method,
      .options = "v"
    )
    
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    .log_msg(sprintf("  Rank %d completed in %.1f mins", k, elapsed), "INFO")
    
    results[[as.character(k)]] <- nmf_fit
    cophenetic_vals[i] <- nmf_fit$cophenetic
    dispersion_vals[i] <- nmf_fit$dispersion
  }
  
  cophenetic_diff <- abs(diff(cophenetic_vals))
  optimal_idx <- which.max(cophenetic_diff) + 1
  optimal_rank <- rank_range[optimal_idx]
  
  .log_msg(sprintf("Optimal rank: k = %d (cophenetic drop = %.3f)", 
                   optimal_rank, cophenetic_diff[optimal_idx - 1]), "INFO")
  
  rank_results <- list(
    fits = results,
    ranks = rank_range,
    cophenetic = cophenetic_vals,
    dispersion = dispersion_vals,
    optimal_rank = optimal_rank,
    cophenetic_diff = cophenetic_diff
  )
  
  dir.create("cache/NMF", showWarnings = FALSE, recursive = TRUE)
  save(rank_results, file = cache_file, compress = TRUE)
  
  return(rank_results)
}

# :: Final NMF Run ------------------------------------------------------------

.nmf_final_run <- function(beta_matrix, rank, nrun = 50, seed = 42, 
                           method = "brunet") {
  .log_msg(sprintf("Running final NMF with rank = %d...", rank), "INFO")
  
  start_time <- Sys.time()
  
  nmf_fit <- nmf(
    t(beta_matrix),
    rank = rank,
    nrun = nrun,
    seed = seed,
    method = method,
    .options = "v"
  )
  
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  .log_msg(sprintf("Final NMF completed in %.1f mins", elapsed), "INFO")
  
  return(nmf_fit)
}

# :: Sample Assignment --------------------------------------------------------

.get_sample_assignments <- function(nmf_fit) {
  assignments <- predict(nmf_fit)
  return(assignments)
}

# :: Visualization Functions ---------------------------------------------------

.plot_cophenetic <- function(rank_results, output_path, config) {
  .log_msg("Plotting cophenetic coefficient...", "INFO")
  
  plot_data <- data.frame(
    rank = rank_results$ranks,
    cophenetic = rank_results$cophenetic,
    diff = c(NA, rank_results$cophenetic_diff)
  )
  
  p <- ggplot(plot_data, aes(x = rank, y = cophenetic)) +
    geom_line(color = "#4EA9E5", size = 1.2) +
    geom_point(color = "#A14444", size = 3, fill = "white", shape = 21) +
    geom_vline(xintercept = rank_results$optimal_rank, 
               linetype = "dashed", color = "#FF6B6B", linewidth = 1) +
    labs(
      title = "NMF Rank Selection via Cophenetic Coefficient",
      subtitle = sprintf("Optimal rank: k = %d", rank_results$optimal_rank),
      x = "Rank (k)",
      y = "Cophenetic Coefficient"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      panel.grid.minor = element_blank()
    )
  
  .save_plot(p, output_path, config)
}

.plot_consensus <- function(nmf_fit, output_path, config) {
  .log_msg("Plotting consensus matrix...", "INFO")
  
  consensus_mat <- consensus(nmf_fit)
  
  Cairo::Cairo(
    file = output_path,
    width = config$viz$width * config$viz$dpi,
    height = config$viz$height * config$viz$dpi,
    type = "png",
    bg = "white"
  )
  
  tinyarray::imageplot(
    consensus_mat,
    Colv = NA,
    Rowv = NA,
    col = Cairo::Cairo.palette(config$viz$consensus_palette),
    main = sprintf("Consensus Matrix (k = %d)", nmf_fit@rank)
  )
  
  dev.off()
  .log_msg(sprintf("Saved: %s", basename(output_path)), "INFO")
}

.plot_clusters <- function(assignments, anno = NULL, output_path, config) {
  .log_msg("Plotting cluster assignments...", "INFO")
  
  plot_data <- data.frame(
    sample = names(assignments),
    cluster = factor(assignments)
  )
  
  if (!is.null(anno)) {
    plot_data <- merge(plot_data, anno, by.x = "sample", by.y = "row.names", all.x = TRUE)
  }
  
  p <- ggplot(plot_data, aes(x = cluster, fill = cluster)) +
    geom_bar() +
    scale_fill_brewer(palette = config$viz$sample_palette) +
    labs(
      title = "Sample Cluster Distribution",
      x = "Cluster",
      y = "Count"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "none"
    )
  
  .save_plot(p, output_path, config)
}

.save_plot <- function(p, base_path, config) {
  dir.create(dirname(base_path), showWarnings = FALSE, recursive = TRUE)
  
  for (fmt in config$viz$formats) {
    output_file <- paste0(base_path, ".", fmt)
    
    if (config$viz$use_cairo && fmt %in% c("png", "tiff")) {
      Cairo::Cairo(
        file = output_file,
        width = config$viz$width * config$viz$dpi,
        height = config$viz$height * config$viz$dpi,
        type = fmt,
        bg = "white"
      )
      print(p)
      dev.off()
    } else {
      ggsave(output_file, p, 
             width = config$viz$width, 
             height = config$viz$height,
             units = config$viz$units,
             dpi = config$viz$dpi)
    }
    
    .log_msg(sprintf("Saved: %s", basename(output_file)), "INFO")
  }
}

# :: Save Results --------------------------------------------------------------

.save_nmf_results <- function(nmf_fit, assignments, output_dir, config) {
  .log_msg("Saving NMF results...", "INFO")
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  if (config$output$file_format == "RData") {
    save(nmf_fit, file = file.path(output_dir, "NMF_run.RData"), 
         compress = config$output$compression)
    .log_msg("Saved: NMF_run.RData", "INFO")
  }
  
  assignments_df <- data.frame(
    sample_id = names(assignments),
    cluster = assignments,
    stringsAsFactors = FALSE
  )
  write.csv(assignments_df, file.path(output_dir, "sample_assignments.csv"), 
            row.names = FALSE)
  .log_msg("Saved: sample_assignments.csv", "INFO")
  
  if (config$output$save_consensus_matrix) {
    consensus <- consensus(nmf_fit)
    save(consensus, file = file.path(output_dir, "consensus_matrix.RData"),
         compress = config$output$compression)
    .log_msg("Saved: consensus_matrix.RData", "INFO")
  }
  
  if (config$output$save_coefficients) {
    coef_list <- list(
      W = basis(nmf_fit),
      H = coef(nmf_fit)
    )
    save(coef_list, file = file.path(output_dir, "nmf_coefficients.RData"),
         compress = config$output$compression)
    .log_msg("Saved: nmf_coefficients.RData", "INFO")
  }
  
  summary <- list(
    timestamp = Sys.time(),
    rank = nmf_fit@rank,
    nrun = nmf_fit@nrun,
    method = nmf_fit@method,
    seed = nmf_fit@seed,
    n_samples = ncol(nmf_fit),
    n_features = nrow(nmf_fit),
    cluster_distribution = table(assignments)
  )
  saveRDS(summary, file.path(output_dir, "nmf_summary.rds"))
  .log_msg("Saved: nmf_summary.rds", "INFO")
}
