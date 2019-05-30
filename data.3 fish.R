

fish.names <- read.table("data/CON fish.tax.txt", header = TRUE, sep="\t")
fish.names <- fish.names[, -1]

fish <- read.table("data/COMSAT 2011 fish.txt", header = TRUE, sep = "\t")
fish <- subset(fish, ID %in% env$ID)

rownames(fish) <- fish$ID
fish <- fish[, -1]

rownames(fish.names) <- fish.names$species
fish.names <- fish.names[, -4]

fish.fam <- as.character(fish.names$family)
fam.tab <- table(fish.fam)
(other.fam <- names(fam.tab[fam.tab < 3]))
fish.fam[fish.fam %in% other.fam] <- "Other"
fish.fam <- factor(fish.fam, 
	levels = c("Cottidae", "Cyprinidae", "Percidae", "Salmonidae", "Other"))
summary(fish.fam)

fish.col <- c(rainbow(length(levels(fish.fam)) - 1, 0.7), gray(0.5, 0.7))
