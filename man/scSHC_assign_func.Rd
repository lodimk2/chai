% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sc_clustering_functions.r
\name{scSHC_assign_func}
\alias{scSHC_assign_func}
\title{Run scSHC}
\usage{
scSHC_assign_func(sce, n_cores = 1)
}
\arguments{
\item{sce}{SingleCellExperiment object of the dataset}

\item{n_cores}{CPU cores for parallel library allocated for running functions.}
}
\description{
Function to run scSHC clustering algorithm.
}
\examples{
\dontrun{
scSHC_assign_func(sce)
}
}
\keyword{scSHC}
