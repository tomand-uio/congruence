
# install.packages(c("maps", "mapdata", "sp", "raster"))
# install.packages(c("rgdal", "tripack", "spdep"))

library(maps)
library(mapdata)
library(sp)
library(raster)
library(rgdal)
library(tripack)
library(spdep)

lake <- read.table("data/COMSAT 2011 lakes.txt", header=TRUE, sep="\t")

excluded.lakes <- c("Jølstravatnet", "Halsjøen", "Bergsvannet", "Sperillen")

lake.sp <- SpatialPointsDataFrame(lake[, c("Longitude", "Latitude")], lake, 
	proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

# Triangulation based Gabriel neighborhood
nb.gabriel <- graph2nb(gabrielneigh(coordinates(lake.sp)), sym = TRUE)

# Administrative boundaries (level=0)
nor <- getData('GADM', country='NOR', level=0)
swe <- getData('GADM', country='SWE', level=0)

# Altitude from worldclim
alt.wc <- getData('worldclim', var='alt', res=2.5)

# Crop to S. Scandinavia
scand <- as(extent(4.5, 19.5, 57, 64), 'SpatialPolygons')
crs(scand) <- crs(alt.wc)
alt.scand <- crop(alt.wc, scand)

pdf("Figure 1.pdf")

plot(alt.scand, col=terrain.colors(255, 0.5),
	xlab="Longitude", ylab="Latitude")
plot(swe, add=TRUE, border=gray(0.5, 0.5))
plot(nor, add=TRUE, border=gray(0.5, 0.5))

plot(nb.gabriel, coordinates(lake.sp), points=FALSE, col="#0000FF20", lwd=2, add = TRUE)

points(lake.sp, pch = 19, cex = log10(lake.sp$Area.km2) + 1, col=gray(0.2, 0.8))
points(subset(lake.sp, Lake %in% excluded.lakes), pch=4, col=2)

area.tick <- c(1, 10, 100)
legend(5.5, 58, paste(area.tick,"km2"), pch = 19, col=gray(0.2, 0.8), 
	pt.cex = log10(area.tick) + 1, bty="n")

dev.off()

