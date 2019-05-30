
source("data.0 env.R")
source("data.1 phyto.R")
source("data.2 zoo.R")
source("data.3 fish.R")

S.phyto <- rowSums(phyto > 0)
S.zoo   <- rowSums(zoo > 0)
S.fish  <- rowSums(fish > 0)

B <- data.frame(phyto=S.phyto, zoo=S.zoo, fish=S.fish)
S <- env[, c("Latitude", "Longitude", "Altitude")] # Spatial
E <- env[, c("log.Cond", "log.TOC", "log.TP")]     # Environmental
D <- env[, c("Longitude", "log.TOC", "log.TP")]    # Design
C <- env[, c("Latitude", "Altitude", "log.Cond")]  # Confounding


#####
# 2-way variance partitioning of R2
# in linear models of species richness 
# by local environment (E) and spatial position (S), 
# and both (ES)

lm.varpart <- function(y, E, S) {
	R2 <- matrix(c(
		summary(lm(y ~ ., data = E))$r.squared,
		summary(lm(y ~ ., data = S))$r.squared,
		summary(lm(y ~ ., data = cbind(S, E)))$r.squared
	), ncol=1)

	# Linear equations for variance explained (M V = R2)
	M <- rbind(c(1,1,0), c(0,1,1), c(1,1,1))
	V <- solve(M, R2)
	return(V) 
}

R2.dc <- cbind(lm.varpart(S.phyto, D, C), lm.varpart(S.zoo, D, C), lm.varpart(S.fish, D, C))
R2.dc <- rbind(R2.dc, 1 - colSums(R2.dc)) # Add residual row

rownames(R2.dc) <- c("D", "DC", "C", "resid.")
colnames(R2.dc) <- c("Phyto", "Zoo", "Fish")

#              Phyto         Zoo       Fish
# D       0.48512965  0.45039522 0.37974661
# DC     -0.01083855 -0.01779218 0.05376406
# C       0.04838267  0.19075896 0.13450386
# resid.  0.47732623  0.37663799 0.43198547

# Design variables (D) explain far more of total R2 than constraints (C)
# Hardly any shared variance between them (DC), which makes sense considering
# they are aligned with different PCA axes (cf. Figure S1.5)

R2.lm <- cbind(lm.varpart(S.phyto, E, S), lm.varpart(S.zoo, E, S), lm.varpart(S.fish, E, S))
R2.lm <- rbind(R2.lm, 1 - colSums(R2.lm)) # Add residual row

rownames(R2.lm) <- c("E", "ES", "S", "resid.")
colnames(R2.lm) <- c("Phyto", "Zoo", "Fish")

#             Phyto        Zoo       Fish
# E      0.05304352 0.03676725 0.01826748
# ES     0.28108322 0.28732750 0.27152057
# S      0.18854704 0.29926725 0.27822648
# resid. 0.47732623 0.37663799 0.43198547

# Local environment alone contribution (E) is rather minuscule, even for phytoplankton
# But there is an appreciable fraction confounded with space (ES)


#####
# 2-way variance partitioning of R2 in RDA models of species composition 
# by local environment (E) and spatial position (S), and both (ES)

library(vegan)

(vp.phyto <- varpart(sqrt(phyto), E, S, transfo="hellinger"))
(vp.zoo <- varpart(zoo, E, S, transfo="hellinger"))
(vp.fish <- varpart(fish, E, S, transfo="hellinger"))

R2.vp <- cbind(
	vp.phyto$part$indfract$Adj.R.squared,
	vp.zoo$part$indfract$Adj.R.squared,
	vp.fish$part$indfract$Adj.R.squared)

rownames(R2.vp) <- c("E", "SE", "S", "resid.")
colnames(R2.vp) <- c("Phyto", "Zoo", "Fish")


#####
# 3-way variance partitioning of R2 
# in linear models of species richness 
# by local environment (E) and spatial position (S), 
# and richness on adjacent trophic levels (B)

G <- rbind(
  c(1, 0, 0, 1, 1, 0, 1),  # M1,E__
  c(0, 1, 0, 1, 0, 1, 1),  # M2,_S_
  c(0, 0, 1, 0, 1, 1, 1),  # M3,__B
  c(1, 1, 0, 1, 1, 1, 1),  # M4,ES_
  c(1, 0, 1, 1, 1, 1, 1),  # M5,E_B
  c(0, 1, 1, 1, 1, 1, 1),  # M6,_SB
  c(1, 1, 1, 1, 1, 1, 1))  # M7,ESB


# Fraction of richness variance explained by environment (E), space (S), 
# biotic interactions (B), and all interactions (E+S, E+B, S+B, E+S+B)

M.ESB <- function(y, E, S, B) {
  M.E__ <- lm(y ~ ., data = as.data.frame(E))
  M._S_ <- lm(y ~ ., data = as.data.frame(S))
  M.__B <- lm(y ~ ., data = as.data.frame(B))
  M.ES_ <- lm(y ~ ., data = cbind(E, S))
  M.E_B <- lm(y ~ ., data = cbind(E, B))
  M._SB <- lm(y ~ ., data = cbind(S, B))
  M.ESB <- lm(y ~ ., data = cbind(E, S, B))
  
  return(list(M.E__, M._S_, M.__B, M.ES_, M.E_B, M._SB, M.ESB))
}

# Using only adjacent trophic levels

M.phyto <- M.ESB(S.phyto, E, S, B[, c("zoo")])
R2.phyto <- sapply(M.phyto, function(x) { summary(x)$r.squared })
V.phyto <- solve(G, R2.phyto)
names(V.phyto) <- c("100", "010", "001", "110", "101", "011", "111")

M.zoo <- M.ESB(S.zoo, E, S, B[,c("phyto", "fish")])
R2.zoo <- sapply(M.zoo, function(x) { summary(x)$r.squared })
V.zoo <- solve(G, R2.zoo)
names(V.zoo) <- c("100", "010", "001", "110", "101", "011", "111")

M.fish <- M.ESB(S.fish, E, S, B[, c("zoo")])
R2.fish <- sapply(M.fish, function(x) { summary(x)$r.squared })
V.fish <- solve(G, R2.fish)
names(V.fish) <- c("100", "010", "001", "110", "101", "011", "111")

round(100 *  cbind(V.phyto, V.zoo, V.fish), digits=2)

#     V.phyto V.zoo V.fish
# 100    4.29  2.91   1.76
# 010    9.18 13.73  13.96
# 001    0.27  0.41   0.30
# 110    7.11  2.72   6.63
# 101    1.01  0.76   0.07
# 011    9.68 16.20  13.86
# 111   21.00 26.02  20.52

round(100 * colSums(cbind(V.phyto, V.zoo, V.fish)), 2)
# V.phyto   V.zoo  V.fish 
#   52.54   62.75   57.10 

#####
# Make figure 3

bar.diagram <- function(R2, col=1:3, main="") {
  barplot(R2, col=col, main=main, ylim=c(0,1), ylab="Partial R2")
  
  R2.cum <- apply(rbind(0 * R2[1, ], R2), 2, cumsum)
  R2.mid <- (R2.cum[1:3, ] + R2.cum[2:4, ]) / 2
  
  x <- rep(1.2 * (1:3) - 0.5, each=3)
  y <- as.numeric(R2.mid)
  text(x, y, rep(rownames(R2), 3), col=gray(0, 0.5))
  
  
  R2.tot <- colSums(R2)
  text(1.2 * (1:3) - 0.5, R2.tot + 0.05, cex=0.8,
       paste(round(100 * R2.tot, digits=1), "%"))
}

dc.col <- c(rgb(1,0.5,0,0.3), rgb(1,1,0,0.3), rgb(0.5,1,0,0.3),rgb(1,1,1,0.3))
vp.col <- c(rgb(1,0,0,0.3), rgb(1,1,0,0.3), rgb(0,1,0,0.3),rgb(1,1,1,0.3))

pdf("Figure 3.pdf")

par(mfrow=c(3, 3), pty="s", mar=c(4, 4, 3, 2))

bar.diagram(R2.dc[1:3, ], col=dc.col, main="Richness (D / C)")
text(3.5, 0.95, "A", font=2, cex=1.8)

bar.diagram(R2.lm[1:3, ], col=vp.col, main="Richness (E / S)")
text(3.5, 0.95, "B", font=2, cex=1.8)

bar.diagram(R2.vp[1:3, ], col=vp.col, main="Community (E / S)")
text(3.5, 0.95, "C", font=2, cex=1.8)

source("function 3-way varpart plot.R")

wheel.diagram(V.phyto, "Phyto")
text(0.95, 0.95, "D", font=2, cex=1.8)

wheel.diagram(V.zoo, "Zoo")
text(0.95, 0.95, "E", font=2, cex=1.8)

wheel.diagram(V.fish, "Fish")
text(0.95, 0.95, "F", font=2, cex=1.8)

dev.off()


