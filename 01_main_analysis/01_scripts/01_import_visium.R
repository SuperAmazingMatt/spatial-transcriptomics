# Import a Space Ranger sample directory (the parent of outs/).
import_visium_section <- function(sample) {
  spe <- BiocIO::import(VisiumIO::TENxVisium(
    spacerangerOut = sample$spaceranger_dir[[1]], sample_id = sample$sample_id[[1]],
    processing = "filtered", format = "h5", images = "lowres"))
  stopifnot(inherits(spe, "SpatialExperiment"), nrow(spe) > 0L, ncol(spe) > 0L,
    all(c("ID", "Symbol") %in% names(SummarizedExperiment::rowData(spe))),
    all(c("in_tissue", "array_row", "array_col", "sample_id") %in%
      names(SummarizedExperiment::colData(spe))))
  if (anyNA(spe$in_tissue)) stop("Missing in_tissue labels.")
  if (!all(as.character(spe$sample_id) == sample$sample_id[[1]])) stop("Sample ID mismatch.")
  spe <- spe[, spe$in_tissue == 1]
  rd <- SummarizedExperiment::rowData(spe)
  if (anyNA(rd$ID) || anyDuplicated(rd$ID)) stop("Ensembl IDs must be present and unique.")
  rownames(spe) <- scuttle::uniquifyFeatureNames(rd$ID, rd$Symbol)
  spe$group <- sample$group[[1]]
  spe$timepoint <- sample$timepoint[[1]]
  if (anyDuplicated(colnames(spe))) stop("Spot barcodes must be unique within a section.")
  spe
}
