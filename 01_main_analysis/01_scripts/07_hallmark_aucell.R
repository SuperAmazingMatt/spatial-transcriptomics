# Retrieve human Hallmark sets mapped to mouse, then score a common Ensembl universe.
# AUCell values are normalized rank-based enrichment scores, not percentages.
score_hallmark_sections <- function(sections, settings) {
  if (!length(sections) || is.null(names(sections)) || anyDuplicated(names(sections))) stop("Name every section uniquely.")
  db <- msigdbr::msigdbr(species = "Mus musculus", collection = "H")
  required <- c("gs_name", "ensembl_gene")
  if (!all(required %in% names(db))) stop("Unsupported msigdbr schema.")
  gene_sets <- lapply(split(db$ensembl_gene, db$gs_name), function(g)
    sort(unique(g[!is.na(g) & nzchar(g)])))
  ids <- lapply(sections, function(s) as.character(SummarizedExperiment::rowData(s)$ID))
  if (any(vapply(ids, function(x) anyNA(x) || anyDuplicated(x) > 0L, logical(1)))) stop("Invalid Ensembl IDs.")
  common <- Reduce(intersect, ids)
  if (!length(common)) stop("No common genes.")
  symbols <- lapply(seq_along(sections), function(i)
    as.character(SummarizedExperiment::rowData(sections[[i]])$Symbol[match(common, ids[[i]])]))
  if (!all(vapply(symbols, identical, logical(1), symbols[[1]]))) stop("Inconsistent common gene annotation.")
  mapped <- lapply(gene_sets, intersect, common)
  coverage <- data.frame(pathway = names(gene_sets), source_genes = lengths(gene_sets),
    measured_genes = lengths(mapped), row.names = NULL)
  mapped <- mapped[lengths(mapped) >= settings$min_genes]
  cutoff <- ceiling(length(common) * settings$rank_fraction)
  if (!length(mapped) || any(lengths(mapped) >= cutoff)) stop("Gene sets incompatible with configured rank cutoff.")
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(settings$seed) # One continuous stream across sections in manifest order.
  scores <- vector("list", length(sections)); names(scores) <- names(sections)
  for (i in seq_along(sections)) {
    x <- as.matrix(SummarizedExperiment::assay(sections[[i]], "logcounts")[match(common, ids[[i]]), , drop = FALSE])
    rownames(x) <- common
    if (any(!is.finite(x))) stop("Non-finite AUCell input.")
    ranks <- AUCell::AUCell_buildRankings(x, featureType = "genes", splitByBlocks = FALSE,
      keepZeroesAsNA = FALSE, BPPARAM = BiocParallel::SerialParam(), plotStats = FALSE, verbose = FALSE)
    auc <- AUCell::AUCell_calcAUC(mapped, ranks, nCores = 1L, normAUC = TRUE,
      aucMaxRank = cutoff, verbose = FALSE)
    scores[[i]] <- as.matrix(AUCell::getAUC(auc))
    stopifnot(identical(colnames(scores[[i]]), colnames(sections[[i]])))
  }
  list(scores = scores, coverage = coverage, common_ensembl = common,
    gene_sets = gene_sets, mapped_gene_sets = mapped, auc_max_rank = cutoff,
    source_versions = if ("db_version" %in% names(db)) unique(db$db_version) else NA_character_,
    seed = settings$seed, sample_order = names(sections))
}
