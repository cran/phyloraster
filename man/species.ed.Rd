% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spp.ed.R
\name{species.ed}
\alias{species.ed}
\title{Calculate Evolutionary distinctiveness for each species}
\usage{
species.ed(tree)
}
\arguments{
\item{tree}{phylo. A dated tree.}
}
\value{
data.frame
}
\description{
This function calculates evolutionary distinctiveness according
to the fair-proportion index for each species.
}
\examples{
\donttest{
library(phyloraster)
tree <- ape::read.tree(system.file("extdata", "tree.nex",
package="phyloraster"))
plot(tree)
ed <- species.ed(tree)
ed
}
}
\references{
Isaac, N. J., Turvey, S. T., Collen, B.,
Waterman, C. and Baillie,
J. E. (2007). Mammals on the EDGE: conservation priorities
based on threat
and phylogeny. PLoS ONE 2, e296.
}
\author{
Neander Marcel Heming and Gabriela Alves-Ferreira
}
