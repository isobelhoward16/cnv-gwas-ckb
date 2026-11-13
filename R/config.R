# Shared helpers: locate the repository root, read the configuration, record step counts.
# All scripts are run from the repository root and source this file first.
suppressPackageStartupMessages({ library(yaml); library(data.table) })
options(width = 200, datatable.print.nrows = 200)

repo_root <- function() {
  d <- normalizePath(getwd())
  repeat {
    if (file.exists(file.path(d, "config", "pipeline.yaml"))) return(d)
    p <- dirname(d)
    if (p == d) stop("Run scripts from the repository root (config/pipeline.yaml not found)")
    d <- p
  }
}

load_config <- function() {
  root <- repo_root()
  pf <- file.path(root, "config", "paths.yaml")
  if (!file.exists(pf)) stop("config/paths.yaml not found: copy config/paths.example.yaml and fill in the paths")
  list(root = root,
       pipeline = read_yaml(file.path(root, "config", "pipeline.yaml")),
       paths = read_yaml(pf))
}

# results/01_callset/<version>/<sub>, created on first use
callset_dir <- function(cfg, ...) {
  d <- file.path(cfg$paths$results_root, "01_callset", cfg$pipeline$callset_version, ...)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

# ---- step counts: every step writes its key numbers to logs/step_counts.tsv ----
counts_file <- function(cfg) file.path(callset_dir(cfg, "logs"), "step_counts.tsv")

start_step <- function(cfg, step) {
  f <- counts_file(cfg)
  if (file.exists(f)) {               # re-running a step replaces its earlier counts
    x <- fread(f, sep = "\t", colClasses = "character")
    fwrite(x[step_name != step], f, sep = "\t")
  }
  cat(sprintf("\n== %s  %s\n", step, format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
}

record <- function(cfg, step, metric, value) {
  v <- if (is.numeric(value)) format(value, big.mark = ",", scientific = FALSE) else as.character(value)
  cat(sprintf("   %-62s %s\n", metric, v))
  f <- counts_file(cfg)
  fwrite(data.table(step_name = step, metric = metric, value = as.character(value)),
         f, sep = "\t", append = file.exists(f) && file.size(f) > 0)
}

check_cols <- function(x, cols, what) {
  miss <- setdiff(cols, names(x))
  if (length(miss)) stop(what, " is missing column(s): ", paste(miss, collapse = ", "))
}
