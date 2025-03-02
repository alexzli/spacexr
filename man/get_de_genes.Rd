% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{get_de_genes}
\alias{get_de_genes}
\title{Returns a list of differentially expressed genes}
\usage{
get_de_genes(
  cell_type_info,
  puck,
  cell_types = NULL,
  fc_thresh = 1.25,
  expr_thresh = 0.00015,
  MIN_OBS = 3
)
}
\arguments{
\item{cell_type_info}{cell type information and profiles of each cell, calculated from the scRNA-seq
reference (see \code{\link{get_cell_type_info}})}

\item{puck}{an object of type \linkS4class{SpatialRNA}}

\item{fc_thresh}{minimum \code{log_e} fold change required for a gene.}

\item{expr_thresh}{minimum expression threshold, as normalized expression (proportion out of 1, or counts per 1).}

\item{MIN_OBS}{the minimum number of occurances of each gene in the SpatialRNA object.}
}
\value{
a list of differntially expressed gene names
}
\description{
For each cell type, chooses genes that have a minimum average normalized expression in that cell
type, and whose expression is larger in that cell type than the average of all cell types.
Filters out mitochondrial genes.
}
