% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/findCellNeighbors.R
\name{findCellNeighbors}
\alias{findCellNeighbors}
\title{Find Cell Neighbors}
\usage{
findCellNeighbors(seurat_obj, slideName_v, radius_v = 50)
}
\arguments{
\item{seurat_obj}{seurat object that has an entry for a SpatialImage object (this will be true if made running buildSeurat)}

\item{slideName_v}{name of the SpatialImage object}

\item{radius_v}{size in um of radius within which to search for neighbors}
}
\description{
Use the spatial coordinates to find all neighbors within a certain radius
}
\details{
Something
}
