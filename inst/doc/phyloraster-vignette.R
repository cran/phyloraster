## ----setup, echo = FALSE, include=FALSE---------------------------------------
library(knitr)
knitr::opts_chunk$set(collapse = TRUE,message=FALSE, warning = FALSE,
                      comment = "#>")

## ----eval = FALSE-------------------------------------------------------------
#   install.packages("phyloraster")

## ----eval = FALSE-------------------------------------------------------------
#  devtools::install_github("gabferreira/phyloraster")

## ----warning = FALSE, message = FALSE-----------------------------------------
library(phyloraster)
library(terra)
library(ape)
library(phylobase)


## -----------------------------------------------------------------------------
data <- load.data.rosauer()
head(data$presab)

## -----------------------------------------------------------------------------
data$tree

## ----fig.height = 5, fig.width = 5, fig.align = 'center'----------------------
plot(data$tree, cex = 0.65)

## -----------------------------------------------------------------------------
data <- load.data.rosauer()
r <- df2rast(x = data$presab, 
             CRS = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
class(r)

## ----fig.height = 5, fig.width = 5, fig.align = 'center'----------------------
plot(r)

## -----------------------------------------------------------------------------
shp <- terra::vect(system.file("extdata", "shps_iucn_spps_rosauer.shp", 
                               package = "phyloraster"))

## ----fig.height = 5, fig.width = 5, fig.align = 'center'----------------------
colors <- rainbow(length(unique(shp$BINOMIAL)),
                  alpha = 0.5)
position <- match(shp$BINOMIAL,
                  unique(shp$BINOMIAL))
colors <- colors[position]
plot(shp, col = colors, lty = 0,
     main = "Spatial polygons")
library(maps)
maps::map(add = TRUE)

## ----message = F--------------------------------------------------------------
r2 <- shp2rast(shp, sps.col = "BINOMIAL", ymask = FALSE, background = 0, 
               resolution = 0.5)
r2
plot(r2[[9]])

## ----message = F--------------------------------------------------------------
library(terra)

shp <- terra::vect(system.file("extdata", "shps_iucn_spps_rosauer.shp",
                              package="phyloraster"))

# create a polygon to use as mask with an extent
e <- terra::ext(113, 123, -43.64, -33.90)
p <- terra::as.polygons(e, crs="")
# cut by the total extension of the polygons
coun.crop <- terra::crop(p, 
                         terra::ext(shp)) 
coun.rast <- terra::rasterize(coun.crop,
terra::rast(terra::ext(shp), resolution = 0.5))

# rasterizing with the mask of the polygon
shp.t <- shp2rast(shp, y = coun.rast, sps.col = "BINOMIAL", ymask = TRUE)
plot(shp.t[[1]], col = c("grey", "green"))

## -----------------------------------------------------------------------------
data <- load.data.rosauer()
names(data$raster) == data$tree$tip.label

## -----------------------------------------------------------------------------
ras <- terra::rast(system.file("extdata", "rast.presab.tif", 
                               package = "phyloraster"))
tree <- ape::read.tree(system.file("extdata", "tree.nex", 
                                   package = "phyloraster"))
dataprep <- phylo.pres(x = ras, tree = tree)

## -----------------------------------------------------------------------------
names(dataprep$x) == tree$tip.label

## ----fig.height = 5, fig.width = 4, fig.align = 'center', warning= FALSE, echo = FALSE----
knitr::include_graphics("figs/tree.jpg")

## ----warning= FALSE-----------------------------------------------------------
ras <- terra::rast(system.file("extdata", "rast.presab.tif", 
                               package = "phyloraster"))
sr <- rast.sr(x = ras)
sr

## ----plot, fig.height = 5, fig.width = 7, fig.align = 'center', warning= FALSE----
plot(sr, main = "Species richness")

## ----warning= FALSE-----------------------------------------------------------
ras <- terra::rast(system.file("extdata", "rast.presab.tif", 
                               package = "phyloraster"))
wer <- rast.we(x = ras)
wer

## ----fig.height = 5, fig.width = 7, fig.align = 'center', warning= FALSE------
wer$WE
plot(wer$WE, main ="Weigthed Endemism")

## ----pdr, warning= FALSE------------------------------------------------------
ras <- terra::rast(system.file("extdata", "rast.presab.tif", 
                               package = "phyloraster"))
tree <- ape::read.tree(system.file("extdata", "tree.nex", 
                                   package = "phyloraster"))
dataprep <- phylo.pres(x = ras, tree = tree, pruning = "tree")

pdr <- rast.pd(x = dataprep$x, edge.path = dataprep$edge.path, 
               branch.length = dataprep$branch.length)

## ----pdr-plot, fig.height = 5, fig.width = 7, fig.align = 'center', warning= FALSE----
plot(pdr$PD, main = "Phylogenetic diversity")

## ----warning= FALSE-----------------------------------------------------------
ras <- terra::rast(system.file("extdata", "rast.presab.tif", 
                               package = "phyloraster"))
tree <- ape::read.tree(system.file("extdata", "tree.nex", 
                                   package = "phyloraster"))
per <- rast.pe(x = dataprep$x, tree)
per

## ----per-plot, fig.height = 5, fig.width = 7, fig.align = 'center', warning= FALSE----
plot(per$PE, main = "Phylogenetic Endemism")

## ----warning= FALSE-----------------------------------------------------------
x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
# phylogenetic tree
tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))
data <- phylo.pres(x, tree)
ed <- rast.ed(data$x, tree)
ed

## ----edr-plot, fig.height = 5, fig.width = 7, fig.align = 'center', warning= TRUE----
terra::plot(ed, main = "Evolutionary Distinctiveness")

## ----warning= FALSE-----------------------------------------------------------
tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))

ed <- phyloraster::species.ed(tree)
head(ed)

## -----------------------------------------------------------------------------
library(SESraster)
ras <- terra::rast(system.file("extdata", "rast.presab.tif", 
                               package = "phyloraster"))
tree <- ape::read.tree(system.file("extdata", "tree.nex", 
                                   package = "phyloraster"))
data <- phylo.pres(ras, tree, pruning = "tree")

t <- rast.pd.ses(data$x, edge.path = data$edge.path, 
                 branch.length = data$branch.length, aleats = 10, 
                 random = "spat")

## ----pds-plot, fig.height = 5, fig.width = 7, fig.align = 'center', warning= FALSE----
plot(t)

## ----sr-plot------------------------------------------------------------------
# load the data
x <- terra::rast(system.file("extdata", "rast.presab.tif", 
                             package="phyloraster"))
# richness
riq.pres <- rast.sr(x)
plot(riq.pres)

## ----srf-plot-----------------------------------------------------------------
# load the data
x <- terra::rast(system.file("extdata", "rast.presab.tif", 
                             package="phyloraster"))
# richness future
riq.fut <- rast.sr(x[[c(1:15)]]) # imagine we lost some species in the future
terra::plot(riq.fut)

## ----dg-----------------------------------------------------------------------
dg <- delta.grid(riq.pres, riq.fut)
plot(dg)

