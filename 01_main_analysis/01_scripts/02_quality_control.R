# Combine global thresholds and the original local-outlier calls.
quality_control_section <- function(spe, settings) {
  mito <- grepl("^mt-", rownames(spe))
  spe <- scuttle::addPerCellQCMetrics(spe, subsets = list(Mito = mito))
  global <- spe$sum < settings$min_umis | spe$detected < settings$min_detected |
    spe$subsets_Mito_percent > settings$max_mito_percent
  spe <- SpotSweeper::localOutliers(spe, metric = "sum", direction = "lower", log = TRUE)
  spe <- SpotSweeper::localOutliers(spe, metric = "detected", direction = "lower", log = TRUE)
  spe <- SpotSweeper::localOutliers(spe, metric = "subsets_Mito_percent",
    direction = "higher", log = FALSE)
  local <- as.logical(spe$sum_outliers) | as.logical(spe$detected_outliers) |
    as.logical(spe$subsets_Mito_percent_outliers)
  if (anyNA(global) || anyNA(local)) stop("Missing QC flags require review.")
  flags <- data.frame(spot_id = colnames(spe), global_outlier = global,
    local_outlier = local, discard = global | local)
  spe <- spe[, !flags$discard]
  if (!ncol(spe)) stop("QC removed every spot.")
  nonzero <- Matrix::rowSums(SummarizedExperiment::assay(spe, "counts")) > 0
  spe <- spe[nonzero, ]
  if (!nrow(spe)) stop("No expressed genes remain.")
  list(spe = spe, flags = flags,
    local_outlier_defaults = formals(SpotSweeper::localOutliers))
}
