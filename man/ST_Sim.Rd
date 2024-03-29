% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chai_st.r
\name{ST_Sim}
\alias{ST_Sim}
\title{Creates a binary similarity matrix from spatial transcriptomics coordinates based on GraphST's method}
\usage{
ST_Sim(spat_data)
}
\arguments{
\item{spat_data}{NAME: the cell  - X: x coordinates  - Y: y coordinates - Cluster: cell type annotation. Ensure that "X" coordinates and "Y" coordinates columns are named with "X" and "Y".}
}
\value{
adj
}
\description{
Creates a binary similarity matrix from spatial transcriptomics coordinates based on GraphST's method
}
