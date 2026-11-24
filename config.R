# Portable method parameters. Sourcing this file does not run an analysis.
project_config <- function(root = getwd()) {
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  list(
    root = root,
    results_dir = Sys.getenv("SPATIAL_RESULTS_DIR", unset = file.path(root,
      "01_main_analysis", "02_results", "generated")),
    qc = list(min_umis = 400, min_detected = 300, max_mito_percent = 20),
    hvg_n = 1000L,
    banksy = list(k_geom = 20L, lambda = 0.8, npcs = 10L, seed = 112358L),
    cluster = list(k = 20L, resolution = 0.6, seed = 123L),
    spanorm = list(seed = 123L, library_target = 10000,
      gene.model = "nb", adj.method = "logpac", sample.p = 0.25,
      scale.factor = 1, df.tps = 6, lambda.a = 1e-4,
      batch = NULL, tol = 1e-4, step.factor = 0.5,
      maxit.nb = 50L, maxit.psi = 25L, maxn.psi = 500L),
    programmes = list(Naive = c("Ccr7", "Sell", "Tcf7"),
      Effector = c("Gzma", "Gzmb", "Bhlhe40"),
      Memory = c("Il7r", "Cd27", "Tcf7")),
    aucell = list(seed = 123L, rank_fraction = 0.05, min_genes = 5L,
      selected = c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_DNA_REPAIR")),
    svg = list(seed = 20261106L, min_count = 3L, min_fraction = 0.005,
      n_neighbors = 10L, order = "AMMD", fdr = 0.05)
  )
}

require_method_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install compatible packages first: ", paste(missing, collapse = ", "))
  invisible(TRUE)
}

validate_sample_manifest <- function(path, root) {
  m <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("sample_id", "spaceranger_dir", "group", "timepoint")
  if (!all(required %in% names(m)) || !nrow(m)) stop("Invalid sample manifest.")
  if (anyNA(m[, required]) || any(!nzchar(as.matrix(m[, required])))) stop("Empty manifest entries.")
  if (anyDuplicated(m$sample_id) || any(!grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", m$sample_id))) {
    stop("Use unique sample IDs containing letters, digits, dot, underscore or hyphen.")
  }
  if (any(grepl("PATH_TO", m$spaceranger_dir, fixed = TRUE))) stop("Replace template paths before running.")
  is_absolute <- grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)", m$spaceranger_dir)
  m$spaceranger_dir[!is_absolute] <- file.path(root, m$spaceranger_dir[!is_absolute])
  m$spaceranger_dir <- vapply(m$spaceranger_dir, normalizePath, character(1),
    winslash = "/", mustWork = TRUE)
  m
}
