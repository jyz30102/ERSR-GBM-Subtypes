#' =============================================================================
#' Module: utils.R
#' Purpose: Utility functions for NMF analysis
#' =============================================================================

.load_config <- function(config_path = "NMF/config.yaml") {
  if (!file.exists(config_path)) {
    stop(sprintf("Config file not found: %s", config_path))
  }
  config <- yaml::read_yaml(config_path)
  return(config)
}

.log_msg <- function(msg, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] [%s] %s\n", timestamp, level, msg))
}

.load_beta_matrix <- function(input_path) {
  .log_msg(sprintf("Loading data: %s", input_path), "INFO")
  
  if (!file.exists(input_path)) {
    stop(sprintf("Input file not found: %s", input_path))
  }
  
  env <- new.env()
  load(input_path, envir = env)
  
  if (exists("betas", envir = env)) {
    beta_matrix <- env$betas
  } else if (exists("GBMbeta", envir = env)) {
    beta_matrix <- env$GBMbeta
  } else {
    stop("No beta matrix found in input file")
  }
  
  .log_msg(sprintf("Loaded: %d samples × %d probes", 
                   nrow(beta_matrix), ncol(beta_matrix)), "INFO")
  
  return(beta_matrix)
}

.validate_beta_matrix <- function(beta_matrix) {
  if (any(beta_matrix < 0, na.rm = TRUE)) {
    warning("Beta matrix contains negative values. Clipping to 0.")
    beta_matrix[beta_matrix < 0] <- 0
  }
  
  if (any(is.na(beta_matrix))) {
    warning("Beta matrix contains NA values. Imputing with column median.")
    for (j in seq_len(ncol(beta_matrix))) {
      na_idx <- which(is.na(beta_matrix[, j]))
      if (length(na_idx) > 0) {
        beta_matrix[na_idx, j] <- median(beta_matrix[, j], na.rm = TRUE)
      }
    }
  }
  
  return(beta_matrix)
}
