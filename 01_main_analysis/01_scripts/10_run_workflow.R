# Sourcing this file only defines a function. Call run_main() explicitly.
# This template has been checked statically, not executed on research data.
run_main <- function(manifest, root = getwd()) {
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  local_environment <- environment()
  sys.source(file.path(root, "config.R"), envir = local_environment)
  script_dir <- file.path(root, "01_main_analysis", "01_scripts")
  for (script in c("01_import_visium.R", "02_quality_control.R", "03_standard_normalisation.R",
      "04_spatial_clustering.R", "05_spanorm.R", "06_marker_programmes.R", "07_hallmark_aucell.R",
      "08_targeted_nnsvg.R", "09_cluster_profiles.R")) {
    sys.source(file.path(script_dir, script), envir = local_environment)
  }
  cfg <- project_config(root)
  samples <- validate_sample_manifest(manifest, root)
  packages <- c("BiocIO", "VisiumIO", "SpatialExperiment", "SummarizedExperiment",
    "SingleCellExperiment", "scuttle", "scran", "SpotSweeper", "Banksy", "igraph",
    "SpaNorm", "Matrix", "msigdbr", "AUCell", "BiocParallel", "nnSVG")
  require_method_packages(packages)
  # Refuse to overwrite an earlier run. Select a fresh output directory.
  if (dir.exists(cfg$results_dir)) stop("Output directory exists; set SPATIAL_RESULTS_DIR to a new directory.")
  dir.create(cfg$results_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(cfg$results_dir)) stop("Cannot create output directory.")
  write_table <- function(x, filename) utils::write.csv(x,
    file.path(cfg$results_dir, filename), row.names = FALSE, na = "")
  saveRDS(cfg, file.path(cfg$results_dir, "parameters.rds"))
  write_table(samples, "sample_manifest.csv")
  write_table(data.frame(package = packages, version = vapply(packages,
    function(p) as.character(utils::packageVersion(p)), character(1))), "package_versions.csv")
  writeLines(capture.output(sessionInfo()), file.path(cfg$results_dir, "session_info.txt"))
  on.exit(writeLines(capture.output(sessionInfo()), file.path(cfg$results_dir, "session_info.txt")), add = TRUE)
  clustered <- adjusted <- programmes <- setNames(vector("list", nrow(samples)), samples$sample_id)
  for (i in seq_len(nrow(samples))) {
    sample_id <- samples$sample_id[i]
    imported <- import_visium_section(samples[i, , drop = FALSE])
    qc <- quality_control_section(imported, cfg$qc)
    write_table(qc$flags, paste0(sample_id, "_qc_flags.csv"))
    writeLines(capture.output(qc$local_outlier_defaults),
      file.path(cfg$results_dir, paste0(sample_id, "_local_outlier_defaults.txt")))
    standard <- standard_normalise_section(qc$spe, cfg$hvg_n)
    clustering <- cluster_section(standard$spe, cfg$banksy, cfg$cluster)
    clustered[[sample_id]] <- clustering$spe
    saveRDS(list(spe = clustering$spe, graph = clustering$graph, variance = standard$variance,
      hvg = standard$hvg), file.path(cfg$results_dir, paste0(sample_id, "_standard_clusters.rds")))
    adjusted[[sample_id]] <- fit_spanorm_section(qc$spe, cfg$spanorm)
    programmes[[sample_id]] <- programme_expression(adjusted[[sample_id]], cfg$programmes)
    saveRDS(adjusted[[sample_id]], file.path(cfg$results_dir, paste0(sample_id, "_spanorm.rds")))
    write_table(data.frame(spot_id = rownames(programmes[[sample_id]]),
      programmes[[sample_id]], check.names = FALSE), paste0(sample_id, "_programmes.csv"))
  }
  hallmark <- score_hallmark_sections(adjusted, cfg$aucell)
  saveRDS(hallmark, file.path(cfg$results_dir, "hallmark_aucell.rds"))
  write_table(hallmark$coverage, "hallmark_gene_coverage.csv")
  svg <- test_programme_genes(adjusted, cfg$programmes, cfg$svg)
  write_table(svg, "programme_gene_nnsvg.csv")
  profiles <- summarise_cluster_profiles(clustered, programmes, hallmark$scores, cfg$aucell$selected)
  write_table(profiles$cluster_means, "cluster_means.csv")
  write_table(profiles$correlations, "cluster_spearman_descriptive.csv")
  writeLines("Workflow completed. Outputs require scientific and diagnostic review before interpretation.",
    file.path(cfg$results_dir, "RUN_STATUS.txt"))
  invisible(list(config = cfg, samples = samples, clustered = clustered, adjusted = adjusted,
    programmes = programmes, hallmark = hallmark, svg = svg, profiles = profiles))
}
