#' Calculates weight matrix from doublet mode results
#'
#' Used during \code{\link{initialize.subtypes}} when unsupervised was run on doublet mode
#'
#' @param RCTD a \code{\linkS4class{RCTD}} object
#' @return a matrix of doublet mode weights
#' @export
weights_from_results <- function(RCTD) {
  df <- RCTD@results$results_df
  wd <- RCTD@results$weights_doublet
  singlets <- rownames(df[df$spot_class == 'singlet', ])
  doublets <- rownames(df[df$spot_class == 'doublet_certain', ])
  barcodes <- c(singlets, doublets)
  cell_types <- levels(df$first_type)
  weights <- matrix(0, nrow = length(barcodes), ncol = length(cell_types))
  rownames(weights) <- barcodes; colnames(weights) <- cell_types
  for (b in singlets) {
    weights[b, as.character(df[b, 'first_type'])] <- 1
  }
  for (b in doublets) {
    weights[b, as.character(df[b, 'first_type'])] <- wd[b, 'first_type']
    weights[b, as.character(df[b, 'second_type'])] <- wd[b, 'second_type']
  }
  return(weights)
}

#' Calculates the weight change between two \code{\linkS4class{RCTD}} objects
#'
#' Used during \code{\link{iter.optim}} as the termination criterion
#'
#' @param RCTD1 a \code{\linkS4class{RCTD}} object
#' @param RCTD2 a \code{\linkS4class{RCTD}} object
#' @return \code{cell_type_info}, which is ready to run the \code{\link{create.RCTD.unsupervised}} function
#' @export
weights_change <- function(RCTD1, RCTD2) {
  weights1 <- as.matrix(RCTD1@results$weights / rowSums(RCTD1@results$weights))
  weights2 <- as.matrix(RCTD2@results$weights / rowSums(RCTD2@results$weights))
  return(norm(weights1 - weights2) / dim(weights1)[1])
}

#' Generates \code{cell_type_info} from cluster assignments on a \code{\linkS4class{SpatialRNA}} object
#'
#' @param puck a \code{\linkS4class{SpatialRNA}} object
#' @param clusters a \code{data.frame} with barcodes as rownames and cluster assignments under a \code{cell_types} column
#' @return \code{cell_type_info}, which is ready to run the \code{\link{create.RCTD.unsupervised}} function
#' @export
cell_type_info_from_clusters <- function(puck, clusters) {
  counts <- puck@counts[, rownames(clusters)]
  nUMI <- puck@nUMI[rownames(clusters)]
  cell_types <- clusters$cell_types
  names(cell_types) <- rownames(clusters)
  cell_type_info <- get_cell_type_info(counts, cell_types, nUMI)
  return(list(info = cell_type_info, renorm = cell_type_info))
}

#' Generates \code{cell_type_info} from CSIDE on a fitted \code{\linkS4class{RCTD}} object
#'
#' @param RCTD a \code{\linkS4class{RCTD}} object to run CSIDE
#' @param cell_types the cell types used for CSIDE. If null, all cell types in \code{cell_type_info} will be chosen.
#' @param cell_type_threshold (Default 10) minimum number of singlets required per cell type. 
#' @param gene_list the genes to compute expression for. If null, all genes will be included.
#' @return \code{cell_type_info}, which is ready to run the \code{\link{create.RCTD.unsupervised}} function
#' @export
cell_type_info_from_de <- function(RCTD, cell_types = NULL, gene_list = NULL, cell_type_threshold = 25) {
  if (is.null(gene_list))
    gene_list <- rownames(RCTD@originalSpatialRNA@counts)
  RCTD@originalSpatialRNA <- restrict_counts(RCTD@originalSpatialRNA, gene_list, UMI_thresh = 0, UMI_max = Inf)
  barcodes <- intersect(names(RCTD@spatialRNA@nUMI), colnames(RCTD@spatialRNA@counts))
  X <- as.matrix(rep(1, length(barcodes)))
  rownames(X) <- barcodes
  if (is.null(cell_types)) {
    cell_type_count <- aggregate_cell_types(RCTD, barcodes, doublet_mode = (RCTD@config$RCTDmode == "doublet"))
    cell_types <- names(which(cell_type_count >= cell_type_threshold))
  }
  RCTD <- run.CSIDE(RCTD, X, barcodes, cell_types, doublet_mode = (RCTD@config$RCTDmode == "doublet"), 
    cell_type_threshold = cell_type_threshold, gene_threshold = -1, sigma_gene = F, test_genes_sig = F, params_to_test = 1)
  cell_type_info <- list(as.data.frame(exp(RCTD@de_results$gene_fits$mean_val)), cell_types, length(cell_types))
  return(list(info = cell_type_info, renorm = cell_type_info))
}

#' Generates \code{cell_type_info} from singlet means on a fitted \code{\linkS4class{RCTD}} object from the original \code{\linkS4class{SpatialRNA}} object
#'
#' @param RCTD a \code{\linkS4class{RCTD}} object
#' @param cell_types the cell types used for CSIDE. If null, all cell types in \code{cell_type_info} will be chosen.
#' @param cell_type_threshold (Default 10) minimum number of singlets required per cell type. 
#' @param gene_list the genes to compute expression for. If null, all genes will be included.
#' @return \code{cell_type_info}, which is ready to run the \code{\link{create.RCTD.unsupervised}} function
#' @export
cell_type_info_from_singlets <- function(RCTD, cell_types = NULL, gene_list = NULL, cell_type_threshold = 25) {
  if (is.null(gene_list))
    gene_list <- rownames(RCTD@originalSpatialRNA@counts)
  results <- RCTD@results$results_df
  singlets <- results[results$spot_class == 'singlet',]
  barcodes <- rownames(singlets)
  cell_type_list <- singlets$first_type
  names(cell_type_list) <- barcodes
  if (is.null(cell_types))
    cell_types <- names(which(table(cell_type_list) > cell_type_threshold))
  cell_type_info <- get_cell_type_info(RCTD@originalSpatialRNA@counts[gene_list, barcodes], cell_type_list, 
    RCTD@originalSpatialRNA@nUMI[barcodes], cell_type_names = cell_types)
  return(list(info = cell_type_info, renorm = cell_type_info))
}
