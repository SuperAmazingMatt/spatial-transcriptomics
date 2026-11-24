# Fixed literature-associated marker panels, not data-selected top genes.
# A spot's programme value is the arithmetic mean of its three gene values.
programme_expression <- function(spe, programmes, assay_name = "logcounts") {
  symbols <- as.character(SummarizedExperiment::rowData(spe)$Symbol)
  genes <- unique(unlist(programmes, use.names = FALSE))
  count <- vapply(genes, function(g) sum(!is.na(symbols) & symbols == g), integer(1))
  if (any(count != 1L)) stop("Each programme gene must map to exactly one measured feature: ",
    paste(names(count)[count != 1L], collapse = ", "))
  x <- as.matrix(SummarizedExperiment::assay(spe, assay_name)[match(genes, symbols), , drop = FALSE])
  rownames(x) <- genes
  if (any(!is.finite(x))) stop("Non-finite programme input.")
  out <- do.call(cbind, lapply(programmes, function(g) colMeans(x[g, , drop = FALSE])))
  rownames(out) <- colnames(spe)
  colnames(out) <- names(programmes)
  out
}
