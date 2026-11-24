# Section-specific BANKSY embedding, Jaccard SNN and Leiden communities.
cluster_section <- function(spe, banksy, cluster) {
  stopifnot("logcounts" %in% SummarizedExperiment::assayNames(spe),
    "hvg" %in% names(SummarizedExperiment::rowData(spe)))
  if (ncol(spe) <= max(banksy$k_geom, cluster$k)) stop("Too few spots for configured neighbours.")
  use <- as.logical(SummarizedExperiment::rowData(spe)$hvg)
  if (sum(use) <= banksy$npcs || ncol(spe) <= banksy$npcs) stop("Too few features/spots for BANKSY PCs.")
  coords <- as.matrix(as.data.frame(SummarizedExperiment::colData(spe)[,
    c("array_row", "array_col"), drop = FALSE]))
  if (any(!is.finite(coords))) stop("Invalid array coordinates.")
  tmp <- spe[use, ]
  SpatialExperiment::spatialCoords(tmp) <- coords
  set.seed(banksy$seed)
  tmp <- Banksy::computeBanksy(tmp, assay_name = "logcounts",
    coord_names = c("array_row", "array_col"), k_geom = banksy$k_geom,
    compute_agf = FALSE)
  args <- list(se = tmp, lambda = banksy$lambda, npcs = banksy$npcs, use_agf = FALSE)
  parameters <- names(formals(Banksy::runBanksyPCA))
  if ("lazy" %in% parameters) args$lazy <- FALSE
  if ("assay_name" %in% parameters) args$assay_name <- "logcounts"
  tmp <- do.call(Banksy::runBanksyPCA, args)
  key <- paste0("PCA_M0_lam", banksy$lambda)
  if (!key %in% SingleCellExperiment::reducedDimNames(tmp)) stop("Unexpected BANKSY reduction names.")
  SingleCellExperiment::reducedDim(spe, "PCA_banksy") <-
    SingleCellExperiment::reducedDim(tmp, key)
  set.seed(cluster$seed)
  graph <- scran::buildSNNGraph(spe, use.dimred = "PCA_banksy", type = "jaccard", k = cluster$k)
  leiden_args <- list(graph = graph, objective_function = "modularity")
  argument <- if ("resolution" %in% names(formals(igraph::cluster_leiden)))
    "resolution" else "resolution_parameter"
  leiden_args[[argument]] <- cluster$resolution
  membership <- igraph::membership(do.call(igraph::cluster_leiden, leiden_args))
  stopifnot(length(membership) == ncol(spe))
  spe$cluster <- factor(membership)
  list(spe = spe, graph = graph)
}
