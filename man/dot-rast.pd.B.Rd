% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rast.pd.ses.R
\name{.rast.pd.B}
\alias{.rast.pd.B}
\title{Calculate phylogenetic diversity for each raster cell}
\usage{
.rast.pd.B(x, edge.path, branch.length, filename = "", ...)
}
\arguments{
\item{x}{SpatRaster. A SpatRaster containing presence-absence data (0 or 1)
for a set of species. The layers (species) will be sorted according to the
tree order. See the phylo.pres function.}

\item{edge.path}{matrix representing the paths through the tree from root
to each tip. See \code{\link{phylo.pres}}}

\item{branch.length}{numeric. A Named numeric vector of branch length for
each species. See \code{\link{phylo.pres}}}

\item{filename}{character. Output filename}

\item{...}{additional arguments passed for terra::app}
}
\value{
SpatRaster
}
\description{
Calculate the sum of the branch length for species present in
raster data.
}
\references{
Faith, D. P. (1992). Conservation evaluation and phylogenetic
diversity. Biological conservation, 61(1), 1-10.
}
\author{
Neander Marcel Heming and Gabriela Alves-Ferreira
}
\keyword{internal}
