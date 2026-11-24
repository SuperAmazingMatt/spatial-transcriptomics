# Optional ordinary Spearman p-values from saved cluster means.
# These nominal diagnostics do not model spatial dependence or test an age effect.
nominal_spearman_diagnostic <- function(cluster_means, pathways) {
  stopifnot(all(c("sample_id", "cluster", "Effector", pathways) %in% names(cluster_means)))
  output <- list()
  for (sample_id in unique(cluster_means$sample_id)) for (pathway in pathways) {
    d <- cluster_means[cluster_means$sample_id == sample_id, , drop = FALSE]
    if (anyDuplicated(d$cluster)) stop("Each cluster must appear once per section.")
    a <- d$Effector; b <- d[[pathway]]
    if (length(a) < 3L || any(!is.finite(a)) || any(!is.finite(b)) ||
        stats::sd(a) == 0 || stats::sd(b) == 0) stop("Invalid correlation input.")
    warnings <- character()
    test <- withCallingHandlers(stats::cor.test(a, b, method = "spearman",
      alternative = "two.sided", exact = TRUE), warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
      })
    output[[length(output) + 1L]] <- data.frame(sample_id = sample_id, pathway = pathway,
      n_clusters = length(a), rho = unname(test$estimate), p_nominal = test$p.value,
      ties = anyDuplicated(a) > 0L || anyDuplicated(b) > 0L,
      warnings = paste(unique(warnings), collapse = " | "))
  }
  do.call(rbind, output)
}
