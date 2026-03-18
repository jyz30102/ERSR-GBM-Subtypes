#' =============================================================================
#' Module: utils.R
#' Purpose: Utility functions for the preprocessing pipeline
#' =============================================================================

#' Log message with timestamp
#' @param msg Message to log
#' @param level Log level (INFO, WARN, ERROR)
#' @param config Configuration list (for .log function)
.log_msg <- function(msg, level = "INFO", config = NULL) {
  if (!is.null(config) && !is.null(config$.log)) {
    config$.log(msg, level)
  } else {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cat(sprintf("[%s] [%s] %s\n", timestamp, level, msg))
  }
}

#' Load filter probe lists from files
#' @param filter_dir Directory containing filter files
#' @param filter_config List of filter file names
#' @return Named list of probe vectors to remove
.load_filter_probes <- function(filter_dir, filter_config) {
  filters <- list()
  
  # Load named filter files
  for (name in names(filter_config)) {
    if (!name %in% c("remove_rs_probes", "remove_ch_probes")) {
      file_path <- file.path(filter_dir, filter_config[[name]])
      if (file.exists(file_path)) {
        probes <- read.table(file_path, header = FALSE, stringsAsFactors = FALSE)[[1]]
        filters[[name]] <- probes
        .log_msg(sprintf("Loaded %s: %d probes", name, length(probes)))
      } else {
        warning(sprintf("Filter file not found: %s", file_path))
      }
    }
  }
  
  return(filters)
}

#' Get probes to remove based on filter criteria
#' @param mset_names Row names of MSet object
#' @param filters List of filter probe vectors
#' @param filter_params Parameters for dynamic filters (rs, ch)
#' @return Vector of row indices to remove
.get_remove_indices <- function(mset_names, filters, filter_params) {
  remove_indices <- integer(0)
  
  # Match named filters
  for (name in names(filters)) {
    if (length(filters[[name]]) > 0) {
      idx <- match(filters[[name]], mset_names, nomatch = 0)
      idx <- idx[idx > 0]  # Remove non-matches
      remove_indices <- c(remove_indices, idx)
      .log_msg(sprintf("  %s: %d probes matched", name, length(idx)), "DEBUG")
    }
  }
  
  # Dynamic filters: rs and ch probes
  if (filter_params$remove_rs_probes) {
    rs_idx <- grep("rs", mset_names)
    remove_indices <- c(remove_indices, rs_idx)
    .log_msg(sprintf("  rs probes: %d matched", length(rs_idx)), "DEBUG")
  }
  
  if (filter_params$remove_ch_probes) {
    ch_idx <- grep("ch", mset_names)
    remove_indices <- c(remove_indices, ch_idx)
    .log_msg(sprintf("  ch probes: %d matched", length(ch_idx)), "DEBUG")
  }
  
  # Return unique indices
  unique(remove_indices)
}

#' Calculate Illumina-style beta values
#' @param methy Methylated signal matrix (probes × samples)
#' @param unmethy Unmethylated signal matrix (probes × samples)
#' @param offset Pseudocount to avoid division by zero
#' @return Beta value matrix
.calc_beta_values <- function(methy, unmethy, offset = 100) {
  beta <- methy / (methy + unmethy + offset)
  # Ensure values are in [0, 1]
  beta <- pmin(pmax(beta, 0), 1)
  return(beta)
}

#' Save R object with optional compression
#' @param obj R object to save
#' @param file_path Output file path
#' @param compress Enable compression
.save_rdata <- function(obj, file_path, compress = TRUE) {
  dir.create(dirname(file_path), showWarnings = FALSE, recursive = TRUE)
  save(obj, file = file_path, compress = compress)
  .log_msg(sprintf("Saved: %s (%.2f MB)", 
                   basename(file_path), 
                   file.info(file_path)$size / 1e6))
}

#' Save beta matrix as CSV (optional)
#' @param beta Beta matrix (samples × probes)
#' @param file_path Output CSV path
#' @param anno Optional annotation to prepend
.save_beta_csv <- function(beta, file_path, anno = NULL) {
  dir.create(dirname(file_path), showWarnings = FALSE, recursive = TRUE)
  
  # Transpose to probes × samples for standard format
  beta_t <- t(beta)
  
  if (!is.null(anno) && nrow(anno) == ncol(beta_t)) {
    # Add annotation columns
    output_df <- cbind(anno[rownames(beta_t), , drop = FALSE], as.data.frame(beta_t))
  } else {
    output_df <- as.data.frame(beta_t)
  }
  
  # Add probe ID column
  output_df <- cbind(probe_id = rownames(output_df), output_df)
  
  write.csv(output_df, file_path, row.names = FALSE)
  .log_msg(sprintf("Saved CSV: %s (%.2f MB)", 
                   basename(file_path), 
                   file.info(file_path)$size / 1e6))
}
