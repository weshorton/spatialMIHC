% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summarizeClusterNeighbors.R
\name{summarizeClusterNeighbors}
\alias{summarizeClusterNeighbors}
\title{Summarize Cluster Neighbors}
\usage{
summarizeClusterNeighbors(
  seurat_obj,
  clusters_v,
  clusterCol_v,
  idCol_v = "ObjectNumber",
  nn_matrix
)
}
\arguments{
\item{seurat_obj}{seurat object made by buildSeurat}

\item{clusters_v}{clusters to iterate over. TO DO - get these from seurat object instead of having to supply them}

\item{clusterCol_v}{name of column in seurat metadata to get clusters from}

\item{idCol_v}{cell id column}

\item{nn_matrix}{nn matrix made earlier}
}
\value{
summary table of neighbors within each cluster
}
\description{
Summarize the neighbor composition of each cluster
}
\details{
to do
}
