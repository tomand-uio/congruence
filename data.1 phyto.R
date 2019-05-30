
pp <- read.table("data/COMSAT 2011 phytoplankton.txt", header = TRUE, sep = "\t")
pp <- subset(pp, Station %in% env$ID)

phyto <- tapply(pp$bio_volume..mm3.l.1., 
	list(as.factor(pp$Station), pp$rubin_kode), sum)
phyto[is.na(phyto)] <- 0 # NAs are true zeros

# Include only taxa present in >3 lakes
phyto <- phyto[, colSums(phyto > 0) > 3] # 73 lakes, 208 taxa
phyto <- phyto[as.character(env$ID), ] # Align with env

# Normalize to relative abundance
phyto <- sweep(phyto, 1, rowSums(phyto), "/")

phyto.tax <- read.table("data/CON phyto.tax.txt", header=TRUE)
phyto.tax <- phyto.tax[colnames(phyto), ]

phyto.class <- as.character(phyto.tax$Class)
class.tab <- table(phyto.class)
other.class <- c("", names(class.tab)[class.tab <= 5])
phyto.class[phyto.class %in% other.class] <- "Other"
phyto.class <- factor(phyto.class, 
	levels = c("Bacillariophyceae", "Chlorophyceae", "Chrysophyceae",    
	"Cryptophyceae", "Cyanophyceae", "Dinophyceae", "Synurophyceae",
	"Trebouxiophyceae", "Zygnematophyceae", "Other"))

phyto.col <- c(rainbow(length(levels(phyto.class)) - 1, 0.7), gray(0.5, 0.7))

