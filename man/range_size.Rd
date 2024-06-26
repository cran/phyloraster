% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/range_size.R
\name{range_size}
\alias{range_size}
\title{Calculate range size for a set of species using a raster as input}
\usage{
range_size(x, cellSz, unit = "m", ...)
}
\arguments{
\item{x}{SpatRaster. A SpatRaster containing presence-absence data (0 or 1)
for a set of species. The layers (species) will be sorted according to the
tree order. See the phylo.pres function.}

\item{cellSz}{SpatRaster. A SpatRaster containing cellSize values.
See \code{\link[terra]{cellSize}}}

\item{unit}{character. One of "m", "km", or "ha"}

\item{...}{additional arguments to be passed passed down from a calling
function.}
}
\value{
vector
}
\description{
This function calculate range size in square meters (by default)
for all cells that are not NA. The size of the cells is constant in degrees
but not in square meters, which was considered in the method applied to
calculate the area.
}
\examples{
\donttest{
x <- terra::rast(system.file("extdata", "rast.presab.tif",
package="phyloraster"))
range_size(x[[1:2]], cellSz <- terra::cellSize(x))
}
}
\author{
Gabriela Alves Ferreira and Neander Marcel Heming
}
