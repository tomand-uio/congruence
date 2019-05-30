
lake <- read.table("data/COMSAT 2011 lakes.txt", header=TRUE, sep="\t")
field <- read.table("data/COMSAT 2011 field.data.txt", header=TRUE)
ctd <- read.table("data/COMSAT 2011 CTD.txt", header=TRUE, sep="\t")
chem <- read.table("data/COMSAT 2011 chemistry.txt", header=TRUE)

lake$Date <- as.Date(lake$Date, "%d.%m.%Y")

env <- lake[, c("ID", "Lake", "Date", "Latitude", "Longitude", "Altitude")]
env <- merge(env, field[, c("ID", "Station.depth", "Secchi.depth")])
env <- merge(env, ctd[, c("ID", "Temp")])
env <- merge(env, chem[, c( "ID", "Cond", "pH", "TOC", "TP")])

# subset(env, (Temp < 10) | (pH < 6) | (pH > 8))
#  404 Jølstravatnet 2011-07-23 61.55793   6.40045
#  482   Bergsvannet 2011-08-09 59.58840  10.05265
# 2888      Halsjøen 2011-07-30 60.86397  12.31110

env <- subset(env, (Temp > 10) & (pH > 6) & (pH < 8))

env[, c(8, 10, 12:13)] <- log(env[, c(8, 10, 12:13)])
names(env) <- c("ID", "Lake", "Date", "Latitude", "Longitude", "Altitude", "Depth",
	"log.Secchi", "Temp", "log.Cond", "pH", "log.TOC", "log.TP")
