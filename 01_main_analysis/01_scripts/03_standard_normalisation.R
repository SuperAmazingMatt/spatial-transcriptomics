# This branch supplies standard logcounts/HVGs for spatial clustering.
# SpaNorm expression is computed separately and does not replace this input.
standard_normalise_section <- function(spe, hvg_n = 1000L) {
  spe <- scuttle::logNormCounts(spe)
  dec <- scran::modelGeneVarByPoisson(spe)
  hvg <- scran::getTopHVGs(dec, n = hvg_n)
  SummarizedExperiment::rowData(spe)$hvg <- rownames(spe) %in% hvg
  list(spe = spe, variance = dec, hvg = hvg)
}
