# Describe existing spatial clusters using separate SpaNorm/AUCell annotations.
summarise_cluster_profiles <- function(clustered, programmes, pathway_scores, selected) {
  if (!setequal(names(clustered), names(programmes)) ||
      !setequal(names(clustered), names(pathway_scores))) stop("Section identities do not match.")
  summaries <- list(); correlations <- list()
  for (sample_id in names(clustered)) {
    spe <- clustered[[sample_id]]
    spots <- colnames(spe); cluster <- as.character(spe$cluster)
    expression <- programmes[[sample_id]]
    scores <- pathway_scores[[sample_id]]
    if (anyNA(cluster) || anyDuplicated(spots) || !setequal(spots, rownames(expression)) ||
        !setequal(spots, colnames(scores))) stop("Spot identities do not match.")
    if (!all(selected %in% rownames(scores))) stop("Selected pathways were not scored.")
    values <- cbind(expression[spots, , drop = FALSE], t(scores[selected, spots, drop = FALSE]))
    if (any(!is.finite(values))) stop("Non-finite cluster annotation.")
    averages <- rowsum(values, cluster, reorder = TRUE)
    sizes <- table(cluster)
    averages <- averages / as.numeric(sizes[rownames(averages)])
    summary <- data.frame(sample_id = sample_id, cluster = rownames(averages),
      n_spots = as.integer(sizes[rownames(averages)]), averages,
      row.names = NULL, check.names = FALSE)
    summaries[[sample_id]] <- summary
    for (pathway in selected) {
      a <- summary$Effector; b <- summary[[pathway]]
      valid <- length(a) >= 3L && stats::sd(a) > 0 && stats::sd(b) > 0
      correlations[[length(correlations) + 1L]] <- data.frame(sample_id = sample_id,
        pathway = pathway, n_clusters = length(a), rho = if (valid)
          stats::cor(a, b, method = "spearman") else NA_real_,
        status = if (valid) "descriptive" else "insufficient_variation")
    }
  }
  list(cluster_means = do.call(rbind, summaries),
    correlations = do.call(rbind, correlations))
}
