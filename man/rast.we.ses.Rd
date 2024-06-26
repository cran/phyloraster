% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rast.we.ses.R
\name{rast.we.ses}
\alias{rast.we.ses}
\title{Calculate weighted endemism standardized for species richness}
\usage{
rast.we.ses(
  x,
  inv.R,
  spat_alg = "bootspat_str",
  spat_alg_args = list(rprob = NULL, rich = NULL, fr_prob = NULL),
  aleats = 10,
  filename = "",
  ...
)
}
\arguments{
\item{x}{SpatRaster. A SpatRaster containing presence-absence data (0 or 1)
for a set of species. The layers (species) will be sorted according to the
tree order. See the phylo.pres function.}

\item{inv.R}{SpatRaster. Inverse of range size. See \code{\link{inv.range}}}

\item{spat_alg}{A function with the algorithm implementing the desired
randomization method. It must work with SpatRaster objects. See examples.
Example of functions that work are: \code{\link[SESraster]{bootspat_naive}},
\code{\link[SESraster]{bootspat_str}}, \code{\link[SESraster]{bootspat_ff}}.}

\item{spat_alg_args}{List of arguments passed to the randomization method
chosen in 'spat_alg'. See \code{\link[SESraster]{bootspat_naive}}, \code{\link[SESraster]{bootspat_str}},
\code{\link[SESraster]{bootspat_ff}}}

\item{aleats}{positive integer. A positive integer indicating how many times
the calculation should be repeated.}

\item{filename}{character. Output filename}

\item{...}{additional arguments passed for terra::app}
}
\value{
SpatRaster
}
\description{
Calculates the standardized effect size for weighted endemism.
See Details for more information.
}
\details{
The spatial randomization (spat) keeps the richness exact and
samples
species presences proportionally to their observed frequency (i.e. number
of occupied pixels). The randomization will not assign values to cells with
nodata. The phylogenetic randomization shuffles taxa names across all taxa
included in phylogeny.
}
\examples{
\donttest{
library(terra)
library(SESraster)
x <- terra::rast(system.file("extdata", "rast.presab.tif",
package="phyloraster"))
t <- rast.we.ses(x[[1:10]], aleats = 3)
plot(t)
}
}
\references{
Williams, P.H., Humphries, C.J., Forey, P.L., Humphries, C.J.,
VaneWright, R.I. (1994). Biodiversity, taxonomic relatedness, and endemism
in conservation. In: Systematics and Conservation Evaluation (eds Forey PL,
Humphries CJ, Vane-Wright RI), p. 438. Oxford University Press, Oxford.

Crisp, M., Laffan, S., Linder, H., Monro, A. (2001).
Endemism in the
Australian flora. Journal of Biogeography, 28, 183–198.
}
\seealso{
\code{\link{phylo.pres}},
\code{\link{inv.range}},
\code{\link{geo.phylo.ses}},
\code{\link{rast.ed.ses}},
\code{\link{rast.pd.ses}},
\code{\link{rast.we.ses}},
\code{\link{rast.pe.ses}},
\code{\link[SESraster]{bootspat_str}},
\code{\link[SESraster]{bootspat_naive}},
\code{\link[SESraster]{bootspat_ff}},
\code{\link[SESraster]{SESraster}}
}
\author{
Neander Marcel Heming and Gabriela Alves-Ferreira
}
