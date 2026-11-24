# Historical alternative normalization, retained as a methods template.
# Its values and any earlier tests are distinct from SpaNorm-derived evidence.
# Supply the original total UMI count per spot before restricting gene rows.
cp10k_programme_expression <- function(counts, symbols, total_umi, programmes) {
  stopifnot(length(symbols) == nrow(counts), length(total_umi) == ncol(counts),
    !is.null(colnames(counts)), !is.null(names(total_umi)))
  if (!setequal(names(total_umi), colnames(counts))) stop("Library-size spot IDs do not match.")
  total_umi <- total_umi[colnames(counts)]
  if (any(!is.finite(total_umi) | total_umi <= 0)) stop("Invalid original library sizes.")
  genes <- unique(unlist(programmes, use.names = FALSE))
  indices <- lapply(genes, function(g) which(!is.na(symbols) & symbols == g))
  names(indices) <- genes
  coverage <- data.frame(gene = genes, feature_rows = lengths(indices),
    present = lengths(indices) > 0L)
  x <- matrix(NA_real_, nrow = length(genes), ncol = ncol(counts),
    dimnames = list(genes, colnames(counts)))
  # The historical method sums exact duplicate Symbol rows before normalization.
  for (gene in genes[coverage$present]) {
    raw <- as.matrix(counts[indices[[gene]], , drop = FALSE])
    if (any(!is.finite(raw)) || any(raw < 0)) stop("Invalid counts.")
    x[gene, ] <- colSums(raw)
  }
  x <- log2(1 + sweep(x, 2L, total_umi, "/") * 10000)
  scores <- do.call(cbind, lapply(programmes, function(g) {
    complete <- all(coverage$present[match(g, coverage$gene)])
    if (complete) colMeans(x[g, , drop = FALSE]) else rep(NA_real_, ncol(counts))
  }))
  rownames(scores) <- colnames(counts); colnames(scores) <- names(programmes)
  list(scores = scores, gene_expression = x, coverage = coverage)
}
