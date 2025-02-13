% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RICE.R
\name{run.RICE.subtypes}
\alias{run.RICE.subtypes}
\title{Runs the subtype pipeline on a \code{\linkS4class{RCTD}} object}
\usage{
run.RICE.subtypes(RCTD, n_iter = 50, MIN_CHANGE = 0.001)
}
\arguments{
\item{RCTD}{an \code{\linkS4class{RCTD}} object created using the \code{\link{create.RCTD}} function.}

\item{n_iter}{maximum number of optimization iterations}

\item{MIN_CHANGE}{minimum change required to terminate optimization}
}
\value{
an \code{\linkS4class{RCTD}} object containing the results of the unsupervised algorithm. Please see \code{\linkS4class{RCTD}}
documentation for more information on interpreting the content of the RCTD object.
}
\description{
Equivalent to sequentially running the functions \code{\link{choose_sigma_c}} and \code{\link{iterOptim}} with subtype mode
}
\details{
Fits only pixels belonging to given cell type(s).
}
