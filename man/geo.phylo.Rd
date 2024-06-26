% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geo.phylo.R
\name{geo.phylo}
\alias{geo.phylo}
\title{Calculate phylogenetic community metrics for raster data}
\usage{
geo.phylo(
  x,
  tree,
  inv.R,
  edge.path,
  branch.length,
  n.descen,
  filename = "",
  ...
)
}
\arguments{
\item{x}{SpatRaster. A SpatRaster containing presence-absence data (0 or 1)
for a set of species. The layers (species) will be sorted according to the
tree order. See the phylo.pres function.}

\item{tree}{phylo. A dated tree.}

\item{inv.R}{SpatRaster. Inverse of range size. See \code{\link{inv.range}}}

\item{edge.path}{matrix representing the paths through the tree from root
to each tip. See \code{\link{phylo.pres}}}

\item{branch.length}{numeric. A Named numeric vector of branch length for
each species. See \code{\link{phylo.pres}}}

\item{n.descen}{numeric. A Named numeric vector of number of descendants for
each branch. See \code{\link{phylo.pres}}}

\item{filename}{character. Output filename}

\item{...}{additional arguments passed for terra::app}
}
\value{
SpatRaster with one layer for each metric
}
\description{
Calculate species richness, phylogenetic diversity, evolutionary
distinctiveness,
phylogenetic endemism and weighted endemism using rasters as input.
}
\details{
Community metrics calculated:
\itemize{
\item{Phylogenetic diversity (Faith 1992)}
\item{Species Richness}
\item{Evolutionary distinctiveness by fair-proportion
(Isaac et al. 2007)}
\item{Phylogenetic endemism (Rosauer et al. 2009)}
\item{Weighted endemism (Crisp et al. 2001, Williams et al. 1994)}
}
}
\examples{
\donttest{
library(terra)
library(phyloraster)
x <- terra::rast(system.file("extdata", "rast.presab.tif",
package="phyloraster"))[[1:10]]
tree <- ape::read.tree(system.file("extdata", "tree.nex",
package="phyloraster"))
data <- phylo.pres(x, tree)
inv.R <- inv.range(data$x)
t <- geo.phylo(data$x, inv.R = inv.R, edge.path = data$edge.path,
branch.length = data$branch.length, n.descen = data$n.descendants)
terra::plot(t)
}
}
\references{
Rosauer, D. A. N., Laffan, S. W., Crisp, M. D., Donnellan, S. C.
and Cook, L. G. (2009). Phylogenetic endemism: a new approach for identifying
geographical concentrations of evolutionary history. Molecular ecology,
18(19), 4061-4072.

Faith, D. P. (1992). Conservation evaluation and phylogenetic
diversity. Biological conservation, 61(1), 1-10.

Williams, P.H., Humphries, C.J., Forey, P.L., Humphries, C.J.
and VaneWright, R.I. (1994). Biodiversity, taxonomic relatedness, and
endemism in conservation. In: Systematics and Conservation Evaluation
(eds Forey PL, Humphries C.J., Vane-Wright RI), p. 438. Oxford University
Press, Oxford.

Crisp, M., Laffan, S., Linder, H. and Monro, A. (2001).
Endemism in the Australian flora. Journal of Biogeography, 28, 183–198.

Isaac, N. J., Turvey, S. T., Collen, B., Waterman, C. and
Baillie, J. E. (2007). Mammals on the EDGE: conservation priorities based on
threat and phylogeny. PLoS ONE 2, e296.

Laffan, S. W., Rosauer, D. F., Di Virgilio, G., Miller, J. T.,
González‐Orozco, C. E., Knerr, N., ... & Mishler, B. D. (2016).
Range‐weighted metrics of species and phylogenetic turnover can better
resolve biogeographic transition zones. Methods in Ecology and Evolution,
7(5), 580-588.
}
\seealso{
\code{\link{phylo.pres}}, \code{\link{inv.range}},
\code{\link{rast.ed}}, \code{\link{rast.pd}},
\code{\link{rast.we}}, \code{\link{rast.pe}}, \code{\link{rast.sr}},
\code{\link{geo.phylo.ses}},
}
\author{
Neander Marcel Heming
}
