
zoo.names <- read.table("data/CON zoo.tax.txt", header = TRUE, sep = "\t")

zoo <- read.table("data/COMSAT 2011 zoo.txt", header = TRUE, sep = "\t")
zoo <- subset(zoo, ID %in% env$ID)

lake <- zoo[, 2:7]
rownames(lake) <- zoo$ID

rownames(zoo) <- zoo$ID
zoo <- zoo[, -(1:7)]

# Recalculate to actual counts
zoo <- sweep(zoo, 1, lake$Fraksjon, "*")

# Remove littoral species
pelagic <- !zoo.names$is.littoral

zoo <- zoo[, pelagic]
zoo.names <- zoo.names[pelagic, ]

rownames(zoo.names) <- zoo.names$RUBIN
zoo.names <- zoo.names[, -1]

# Normalize to relative abundance
zoo <- sweep(zoo, 1, rowSums(zoo), "/")

zoo.fam <- as.character(zoo.names$family)
fam.tab <- table(zoo.fam)
(other.fam <- names(fam.tab[fam.tab < 3]))
zoo.fam[zoo.fam %in% other.fam] <- "Other"
zoo.fam <- factor(zoo.fam, 
	levels = c("Bosminidae", "Cyclopidae", "Daphniidae",
		"Diaptomidae", "Temoridae", "Other"))
summary(zoo.fam)

zoo.col <- c(rainbow(length(levels(zoo.fam)) - 1, 0.7), gray(0.5, 0.7))

