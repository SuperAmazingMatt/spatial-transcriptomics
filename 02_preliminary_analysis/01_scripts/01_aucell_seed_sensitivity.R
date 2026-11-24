# Optional check of AUCell tie-breaking, using the same common universe/method.
# Requires score_hallmark_sections() from the main analysis module.
aucell_seed_sensitivity <- function(sections, settings, seeds = c(123L, 124L, 125L)) {
  if (anyDuplicated(seeds)) stop("Use distinct seeds.")
  runs <- lapply(seeds, function(seed) {
    current <- settings; current$seed <- seed
    score_hallmark_sections(sections, current)
  })
  names(runs) <- as.character(seeds)
  # Compare scores only if the retrieved collection and measured universe agree.
  reference <- runs[[1]]
  for (run in runs) stopifnot(identical(run$common_ensembl, reference$common_ensembl),
    identical(run$mapped_gene_sets, reference$mapped_gene_sets))
  comparisons <- list()
  for (i in seq_along(runs)[-1L]) for (sample_id in names(sections)) {
    a <- reference$scores[[sample_id]]; b <- runs[[i]]$scores[[sample_id]]
    stopifnot(identical(dimnames(a), dimnames(b)))
    comparisons[[length(comparisons) + 1L]] <- data.frame(sample_id = sample_id,
      reference_seed = seeds[1], comparison_seed = seeds[i], pathway = rownames(a),
      maximum_absolute_score_change = apply(abs(a - b), 1L, max))
  }
  list(runs = runs, comparisons = if (length(comparisons)) do.call(rbind, comparisons) else data.frame())
}
