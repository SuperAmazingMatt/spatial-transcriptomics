# The resulting logcounts assay already contains log2(1 + PAC).
fit_spanorm_section <- function(spe, settings) {
  counts <- SummarizedExperiment::assay(spe, "counts")
  totals <- as.numeric(Matrix::colSums(counts))
  if (any(!is.finite(totals) | totals <= 0)) stop("Non-positive library size.")
  SingleCellExperiment::sizeFactors(spe) <- totals / settings$library_target
  coords <- SpatialExperiment::spatialCoords(spe)
  ids <- dimnames(spe)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(settings$seed)
  args <- settings[setdiff(names(settings), c("seed", "library_target"))]
  args <- c(list(spe), args, list(overwrite = TRUE, verbose = TRUE))
  fitted <- do.call(SpaNorm::SpaNorm, args)
  stopifnot(identical(dimnames(fitted), ids),
    identical(SpatialExperiment::spatialCoords(fitted), coords),
    isTRUE(all.equal(SummarizedExperiment::assay(fitted, "counts"), counts)),
    isTRUE(all.equal(as.numeric(SingleCellExperiment::sizeFactors(fitted)),
      totals / settings$library_target)))
  x <- SummarizedExperiment::assay(fitted, "logcounts")
  if (any(!is.finite(x)) || any(x < 0)) stop("Invalid SpaNorm adjusted expression.")
  fitted
}
