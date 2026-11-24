# Per-gene spatial tests on SpaNorm expression. The requested family is fixed
# before fitting; sparse/failed requests stay in the multiplicity denominator.
test_programme_genes <- function(sections, programmes, settings) {
  genes <- unique(unlist(programmes, use.names = FALSE))
  family_size <- length(genes) * length(sections)
  if (!family_size || is.null(names(sections))) stop("Named sections and target genes are required.")
  rows <- vector("list", family_size); cursor <- 0L
  for (i in seq_along(sections)) {
    spe <- sections[[i]]
    symbols <- as.character(SummarizedExperiment::rowData(spe)$Symbol)
    counts <- SummarizedExperiment::assay(spe, "counts")
    expression <- SummarizedExperiment::assay(spe, "logcounts")
    coords <- as.matrix(SpatialExperiment::spatialCoords(spe))
    if (nrow(coords) != ncol(spe) || any(!is.finite(coords))) stop("Invalid spatial coordinates.")
    if (ncol(spe) <= settings$n_neighbors) stop("Too few spots for configured nnSVG neighbours.")
    for (gene in genes) {
      cursor <- cursor + 1L
      index <- which(!is.na(symbols) & symbols == gene)
      row <- data.frame(sample_id = names(sections)[i], gene = gene,
        n_spots = ncol(spe), detected_spots = NA_integer_,
        minimum_spots = ceiling(settings$min_fraction * ncol(spe)),
        status = "ready", sigma.sq = NA_real_, tau.sq = NA_real_, phi = NA_real_,
        prop_sv = NA_real_, loglik = NA_real_, loglik_lm = NA_real_, LR_stat = NA_real_,
        p_package = NA_real_, p_stable = NA_real_, log_p_stable = NA_real_,
        seed = settings$seed + i, warnings = "", error = "", stringsAsFactors = FALSE)
      if (length(index) != 1L) {
        row$status <- if (!length(index)) "missing_feature" else "ambiguous_feature"
      } else {
        row$detected_spots <- sum(as.numeric(counts[index, ]) >= settings$min_count)
        x <- as.matrix(expression[index, , drop = FALSE]); rownames(x) <- gene
        if (row$detected_spots < row$minimum_spots) row$status <- "too_sparse"
        else if (any(!is.finite(x))) row$status <- "invalid_expression"
        else if (stats::var(as.numeric(x)) == 0) row$status <- "constant_expression"
      }
      if (row$status == "ready") {
        captured <- character()
        RNGkind("Mersenne-Twister", "Inversion", "Rejection")
        set.seed(row$seed)
        fit <- tryCatch(withCallingHandlers(nnSVG::nnSVG(input = x,
          spatial_coords = coords, X = NULL, n_neighbors = settings$n_neighbors,
          order = settings$order, n_threads = 1L, BPPARAM = BiocParallel::SerialParam(),
          verbose = FALSE), warning = function(w) {
            captured <<- c(captured, conditionMessage(w)); invokeRestart("muffleWarning")
          }), error = function(e) e)
        row$warnings <- paste(unique(captured), collapse = " | ")
        if (inherits(fit, "error")) {
          row$status <- "fit_failed"; row$error <- conditionMessage(fit)
        } else {
          statistics <- c("sigma.sq", "tau.sq", "phi", "prop_sv", "loglik", "loglik_lm", "LR_stat")
          if (nrow(fit) != 1L || !identical(rownames(fit), gene) ||
              !all(c(statistics, "pval") %in% colnames(fit))) stop("Unexpected nnSVG return schema.")
          for (key in statistics) row[[key]] <- as.numeric(fit[1, key])
          row$p_package <- as.numeric(fit[1, "pval"])
          valid <- all(is.finite(unlist(row[, statistics]))) &&
            row$sigma.sq >= 0 && row$tau.sq >= 0 && row$phi > 0 &&
            row$prop_sv >= 0 && row$prop_sv <= 1
          row$status <- if (valid) "tested" else "fit_invalid"
          if (valid) {
            row$p_stable <- stats::pchisq(max(0, row$LR_stat), df = 2, lower.tail = FALSE)
            row$log_p_stable <- -max(0, row$LR_stat) / 2
          }
        }
      }
      rows[[cursor]] <- row
    }
  }
  out <- do.call(rbind, rows)
  tested <- out$status == "tested"
  p <- out$p_stable; p[!tested] <- 1
  out$q_global <- stats::p.adjust(p, method = "BH", n = family_size)
  out$q_global[!tested] <- NA_real_
  # Compute BH in log space as well, retaining very small probabilities.
  lp <- out$log_p_stable; lp[!tested] <- 0
  order_p <- order(lp)
  adjusted <- pmin(0, rev(cummin(rev(lp[order_p] + log(family_size) - log(seq_along(lp))))))
  out$log_q_global <- numeric(family_size); out$log_q_global[order_p] <- adjusted
  out$log_q_global[!tested] <- NA_real_
  out$supported <- tested & !is.na(out$log_q_global) & out$log_q_global < log(settings$fdr)
  out$numeric_underflow <- tested & (out$p_stable == 0 | out$q_global == 0)
  out$variance_boundary <- tested & (out$prop_sv <= 1e-6 | out$prop_sv >= 1 - 1e-6)
  out$evidence <- ifelse(out$status == "too_sparse", "Sparse",
    ifelse(!tested, "Untestable", ifelse(out$supported, "Supported", "No support")))
  out$family_size <- family_size
  out
}
