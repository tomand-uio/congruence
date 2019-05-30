
source("data.0 env.R")
source("data.1 phyto.R")
source("data.2 zoo.R")
source("data.3 fish.R")

S.obs <- data.frame(
	S.phyto = rowSums(phyto > 0),
	S.zoo = rowSums(zoo > 0),
	S.fish = rowSums(fish > 0))

summary(S.obs)

(rho.zoo.phyto <- with(S.obs, cor.test(S.zoo, S.phyto, method="spearman")))   # rho = 0.59, p < 0.001
(rho.fish.phyto <- with(S.obs, cor.test(S.fish, S.phyto, method="spearman"))) # rho = 0.58, p < 0.001
(rho.fish.zoo <- with(S.obs, cor.test(S.fish, S.zoo, method="spearman")))     # rho = 0.70, p < 0.001

library(vegan)

(mds.phyto <- metaMDS(sqrt(phyto), k=2, trymax=1000)) 	# Conv. 38 itn. stress = 0.20
(mds.zoo <- metaMDS(zoo, k=2, trymax=1000))			# Conv. 33 itn. stress = 0.19
(mds.fish <- metaMDS(fish, k=2, trymax=1000))			# Conv.  3 itn. stress = 0.14

(pro.phyto.zoo  <- protest(mds.phyto, mds.zoo))  # corr = 0.51, p = 0.001 on 999 perm.
(pro.phyto.fish <- protest(mds.phyto, mds.fish)) # corr = 0.54, p = 0.001 on 999 perm.
(pro.zoo.fish   <- protest(mds.zoo, mds.fish))   # corr = 0.61, p = 0.001 on 999 perm.

pdf("Figure 2.pdf")

par(mfcol=c(3, 3), pty="s", mar=c(4, 4, 3, 2) + 0.1)

plot(S.zoo ~ S.phyto, pch=19, cex=1.5, data=S.obs, 
	xlim=c(0, 100), xlab="Phytoplankton richness", 
	ylim=c(-2, 25), ylab="Zooplankton richness")
text(50, -1, paste("Spearman corr. =", round(rho.zoo.phyto$estimate , 2)))
text(5, 23, "A", font=2, cex=1.8)

plot(S.fish ~ S.phyto, pch=19, cex=1.5, data=S.obs, 
	xlim=c(0, 100), xlab="Phytoplankton richness", 
	ylim=c(-2, 25), ylab="Fish richness")
text(50, -1, paste("Spearman corr. =", round(rho.fish.phyto$estimate , 2)))
text(5, 23, "B", font=2, cex=1.8)

plot(S.fish ~ S.zoo, pch=19, cex=1.5, data=S.obs, 
	xlim=c(0, 25), xlab="Zooplankton richness", 
	ylim=c(-2, 25), ylab="Fish richness")
text(12.5, -1, paste("Spearman corr. =", round(rho.fish.zoo$estimate , 2)))
text(2, 23, "C", font=2, cex=1.8)

with(pro.phyto.zoo, plot(X[, 1], Yrot[,1], pch=19, cex=1.5,
	xlim=c(-0.3, 0.3), ylim=c(-0.15, 0.15),
	xlab="NMDS1(Phytoplankton)", ylab="Rotated NMDS1(Zooplankton)"))
text(0, -0.14, paste("Procrustes corr. =", round(pro.phyto.zoo$scale, 2)))
text(-0.26, 0.13, "D", font=2, cex=1.8)

with(pro.phyto.fish, plot(X[, 1], Yrot[,1], pch=19, cex=1.5,
	xlim=c(-0.3, 0.3), ylim=c(-0.15, 0.15),
	xlab="NMDS1(Phytoplankton)", ylab="Rotated NMDS1(Fish)"))
text(0, -0.14, paste("Procrustes corr. =", round(pro.phyto.fish$scale, 2)))
text(-0.26, 0.13, "E", font=2, cex=1.8)

with(pro.zoo.fish, plot(X[, 1], Yrot[,1], pch=19, cex=1.5,
	xlim=c(-0.3, 0.3), ylim=c(-0.15, 0.15),
	xlab="NMDS1(Zooplankton)", ylab="Rotated NMDS1(Fish)"))
text(0, -0.14, paste("Procrustes corr. =", round(pro.zoo.fish$scale, 2)))
text(-0.26, 0.13, "F", font=2, cex=1.8)

dev.off()
