#' =============================================================================
#' Module: config.R
#' Purpose: Load and manage project configuration
#' =============================================================================

#' Load configuration from YAML file
#' @param config_path Path to config.yaml (relative to project root)
#' @return List containing configuration settings
.load_config <- function(config_path = "config/config.yaml") {
  # Check if yaml package is available
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' required for configuration. Install with: install.packages('yaml')")
  }
  
  # Resolve path relative to project root
  project_root <- .get_project_root()
  config_file <- file.path(project_root, config_path)
  
  if (!file.exists(config_file)) {
    warning(sprintf("Config file not found: %s. Using defaults.", config_file))
    return(.get_default_config())
  }
  
  config <- yaml::read_yaml(config_file)
  .log(sprintf("Configuration loaded from %s", config_file), "INFO")
  
  # Resolve all paths to absolute
  config$paths <- .resolve_paths(config$paths, project_root)
  
  return(config)
}

#' Get project root directory (directory containing config/)
#' @return Character string of project root path
.get_project_root <- function() {
  # Strategy 1: Check for config/ directory in current or parent dirs
  start_dir <- getwd()
  current <- start_dir
  
  for (i in 1:5) {  # Check up to 5 levels up
    if (dir.exists(file.path(current, "config")) && 
        file.exists(file.path(current, "config", "config.yaml"))) {
      return(current)
    }
    current <- dirname(current)
  }
  
  # Strategy 2: Fallback to script directory
  script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
  .log(sprintf("Using script directory as project root: %s", script_dir), "WARN")
  return(script_dir)
}

#' Resolve relative paths to absolute paths
#' @param paths List of path strings
#' @param base_dir Base directory for resolution
#' @return List with absolute paths
.resolve_paths <- function(paths, base_dir) {
  lapply(paths, function(p) {
    if (is.character(p) && !grepl("^(/|[A-Z]:)", p)) {  # Relative path
      file.path(base_dir, p)
    } else {
      p  # Already absolute
    }
  })
}

#' Get default configuration (fallback)
#' @return List with default settings
.get_default_config <- function() {
  list(
    paths = list(
      raw_idat = "data/raw/GSE90496_RAW",
      filter_dir = "data/filters",
      output_dir = "results/GBMLGG",
      log_dir = "results/GBMLGG/logs",
      functions_script = "R/MNPprocessIDAT_functions.R"
    ),
    geo = list(
      accession = "GSE90496",
      target_classes = c("LGG", "GBM"),
      methylation_class_col = "methylation class:ch1",
      material_col = "material:ch1"
    ),
    preprocessing = list(
      filters = list(
        ambiguous = "amb_3965probes.vh20151030.txt",
        epic = "epicV1B2_32260probes.vh20160325.txt",
        snp = "snp_7998probes.vh20151030.txt",
        xy = "xy_11551probes.vh20151030.txt",
        remove_rs_probes = TRUE,
        remove_ch_probes = TRUE
      ),
      batch_correction = list(
        method = "limma",
        batch_col = "material:ch1",
        batch_levels = c(Frozen = 1, FFPE = 2),
        pseudocount = 1,
        offset = 100
      ),
      beta = list(
        formula = "methy / (methy + unmethy + offset)",
        offset = 100
      )
    ),
    output = list(
      save_intermediate = TRUE,
      save_final_beta = TRUE,
      file_format = "RData",
      compression = TRUE
    ),
    logging = list(
      level = "INFO",
      log_to_file = TRUE,
      log_to_console = TRUE,
      timestamp_format = "%Y-%m-%d %H:%M:%S"
    )
  )
}

#' Initialize logging system
#' @param config Configuration list
#' @return NULL (side effect: sets up logging)
.init_logging <- function(config) {
  log_cfg <- config$logging
  
  # Create log directory
  if (log_cfg$log_to_file) {
    dir.create(config$paths$log_dir, showWarnings = FALSE, recursive = TRUE)
    log_file <- file.path(
      config$paths$log_dir, 
      sprintf("preprocessing_%s.log", format(Sys.time(), "%Y%m%d_%H%M%S"))
    )
    # Simple logging function attached to config
    config$.log <- function(msg, level = "INFO") {
      timestamp <- format(Sys.time(), log_cfg$timestamp_format)
      log_line <- sprintf("[%s] [%s] %s\n", timestamp, level, msg)
      if (log_cfg$log_to_console) cat(log_line)
      if (log_cfg$log_to_file) cat(log_line, file = log_file, append = TRUE)
    }
  } else {
    config$.log <- function(msg, level = "INFO") {
      timestamp <- format(Sys.time(), log_cfg$timestamp_format)
      cat(sprintf("[%s] [%s] %s\n", timestamp, level, msg))
    }
  }
  
  config$.log("Logging initialized", "INFO")
  return(config)
}
