% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RICE_utils.R
\name{weights_from_results}
\alias{weights_from_results}
\title{Calculates weight matrix from doublet mode results}
\usage{
weights_from_results(RCTD)
}
\arguments{
\item{RCTD}{a \code{\linkS4class{RCTD}} object}
}
\value{
a matrix of doublet mode weights
}
\description{
Used during \code{\link{initialize.subtypes}} when unsupervised was run on doublet mode
}
