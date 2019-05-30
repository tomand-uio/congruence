
source("data.0 env.R")
source("data.1 phyto.R")
source("data.2 zoo.R")
source("data.3 fish.R")

env$S.phyto <- rowSums(phyto > 0)
env$S.zoo   <- rowSums(zoo   > 0)
env$S.fish  <- rowSums(fish  > 0)

# Using only design covariates (Longitude,  TOC, total P) since
# they explain most of the species richness variance  (Fig. 3A)

env.sem <- as.data.frame(scale(env[, c(5, 12:16)]))
plot(env.sem)

################################################################################
######        Non-recursive SEM Models for trophic interactions         ######
################################################################################

# Based on script in supplementary data from: 
# Zhang, J., Qian, H., Girardello, M., Pellissier, V., Nielsen, S. E., & Svenning, J. C. (2018). 
# Trophic interactions among vertebrate guilds and plants shape global patterns in species 
# diversity. Proceedings of the Royal Society B: Biological Sciences, 285(1883), 20180949.

# https://figshare.com/articles/Electronic_supplementary_material_R_script_from_Trophic_interactions_among_vertebrate_guilds_and_plants_shape_global_patterns_in_species_diversity/6791282

library(systemfit)

# Assuming that only adjacent trophic levels have effects

equ.list <- list(
  m.phyto = S.phyto ~ Longitude + log.TP + log.TOC + S.zoo,
  m.zoo   = S.zoo   ~ Longitude + log.TP + log.TOC + S.phyto + S.fish,
  m.fish  = S.fish  ~ Longitude + log.TP + log.TOC + S.zoo
)

summary(OLS.fit <- systemfit(equ.list, method="OLS", data=env.sem))

#          N  DF     SSR  detRCov   OLS-R2 McElroy-R2
# system 219 203 108.442 0.128741 0.497955    0.59972

# Significant biotic interactions for zoo <--> fish but not for phyto

# Simultaneous equations fitted with ordinary leasts quares (OLS)
# assume that residuals of individual equations are uncorrelated

#           m.phyto     m.zoo    m.fish
# m.phyto  1.000000 -0.200899  0.189530
# m.zoo   -0.200899  1.000000 -0.322015
# m.fish   0.189530 -0.322015  1.000000

# But since there seems to be certain residual correlations
# it may be be worthwhile to also try a generalized least squares 
# approach (which in systemfit goes under the name "Seemingly
# Unrelated Regressions" (SUR))

summary(SUR.fit <- systemfit(equ.list, method="SUR", data=env.sem))

#          N  DF     SSR detRCov   OLS-R2 McElroy-R2
# system 219 203 115.145  0.0619 0.466922   0.839714

# McElroy's R2 is substantially higher for the SUR/GLS model
# but the assumptions are also more restrictive

# Correcting sample sizes and standard errors in SEMs with 
# spatial autocorrelation in endogenous variables
# https://github.com/jebyrnes/spatial_correction_lavaan

source("function SpatialCorrect.R")

OLS.sp.cor <- lavSpatialCorrect(obj=OLS.fit, xvar=env$Longitude, yvar=env$Latitude)
do.call(rbind, OLS.sp.cor$Morans_I)

SUR.sp.cor <- lavSpatialCorrect(obj=SUR.fit, xvar=env$Longitude, yvar=env$Latitude)
do.call(rbind, SUR.sp.cor$Morans_I)

# Non-significant Moran's I autocorrelation on both model residuals
# means that no corrections in sample size were done and that
# corrected values are identical with original...

#####
# Construct a network from the coefficient table of a systemfit model

coef.network <-  function(model, node.names) {
  
  # Extract coefficients and remove Intercepts term
  
  b <- model$coefficients
  b <- b[!grepl("Intercept",names(b))]
  
  # Construct edge list from coefficient table
  
  b.edge <- sub("m.", "S.", names(b))
  b.edge <- strsplit(b.edge, "_")
  b.edge <- t(matrix(unlist(b.edge), nrow=2))
  
  nodes <- unique(as.character(b.edge))
  to <- factor(b.edge[, 1], levels=nodes) 
  from <- factor(b.edge[, 2], levels=nodes)
  
  # Simplify node labels
  
  to <- factor(node.names[to], levels=node.names) 
  from <- factor(node.names[from], levels=node.names)
  
  return(data.frame(to, from, b))
}
  
nodes <- c("Phyto", "Zoo", "Fish", "Long", "TP", "TOC")
OLS.network <- coef.network(OLS.fit, node.names=nodes) 
SUR.network <- coef.network(SUR.fit, node.names=nodes) 

# Needed for pdf export
# devtools::install_github('rich-iannone/DiagrammeRsvg') 
# install.packages("rsvg")

library(DiagrammeR)
library(DiagrammeRsvg)

# Make DiagrammeR graph object from coefficients and node  lists

sem.graph <- function(to, from, b) {
  edf <- create_edge_df(from=from, to=to, label=round(b, 2), fontsize=7, 
         fontcolor="Gray40", color=c("ForestGreen", "Tomato")[1 + (b < 0)], penwidth=10*abs(b))
  ndf <- create_node_df(n=length(nodes), label=nodes, fontcolor="Gray40",
         fillcolor=rep(c("PowderBlue", "Ivory"), each=3), shape=rep(c("circle", "rectangle"), each=3))
  return(create_graph(nodes_df=ndf, edges_df=edf))
}

OLS.graph <- with(OLS.network, sem.graph(to, from, b))
SUR.graph <- with(SUR.network, sem.graph(to, from, b))

render_graph(OLS.graph)
export_graph(OLS.graph, file_name="Figure 4A.pdf")

render_graph(SUR.graph)
export_graph(SUR.graph, file_name="Figure 4B.pdf")



################################################################################
######        Recursive SEM Models for trophic interactions         ######
################################################################################

# Non-recursive (systemfit ett al.) can have loops or reciprocal effects
# Recursive models (lavaan etc) can only have unidirectional causal effects 

library(lavaan)

# Top-down model:
# Support for Fish -> Zoo, less so for Zoo -> Phyto

M1 <- "
  S.phyto ~ Longitude + log.TP + log.TOC + S.zoo
  S.zoo   ~ Longitude + log.TP + log.TOC + S.fish
  S.fish  ~ Longitude + log.TP + log.TOC
"

summary(sem1 <- sem(M1, data=env.sem)) # P-value (Chi-square) = 0.102 (1 df)
modindices(sem1,  minimum.value=1, free.remove=TRUE) # All <3

lavInspect(sem1, "R2")
# S.phyto   S.zoo  S.fish 
#   0.497   0.493   0.434 

# Bottom-up model:
# Support for Zoo -> Fish, less so for Phyto -> Zoo

M2 <- "
  S.phyto ~ Longitude + log.TP + log.TOC
  S.zoo   ~ Longitude + log.TP + log.TOC + S.phyto
  S.fish  ~ Longitude + log.TP + log.TOC + S.zoo
"

summary(sem2 <- sem(M2, data=env.sem)) # P-value (Chi-square) = 0.102 (1 df)
modindices(sem2,  minimum.value=1, free.remove=TRUE) # All <3

lavInspect(sem2, "R2")
# S.phyto   S.zoo  S.fish 
#   0.474   0.458   0.493 

# Bidirectional zoo-fish interaction model

M3 <- "
  S.phyto ~ Longitude + log.TP + log.TOC
  S.zoo   ~ Longitude + log.TP + log.TOC + S.fish
  S.fish  ~ Longitude + log.TP + log.TOC + S.zoo
"

summary(sem3 <- sem(M3, data=env.sem)) # P-value (Chi-square) = 0.015 (1 df)

# Could not compute standard errors! The information matrix could not be inverted. 
# This may be a symptom that the model is not identified.

# In other words: lavaan cannot fit models with bidirectional zoo - fish interaction
# nor can it identify the "adjacent trophic levels" model fitted by systemfit 

sem.network <-  function(model, node.names) {
  
  # Extract regression coefficients (op == "~")
  
  b <- parTable(model)
  b <- subset(b, op == "~")[,  c("lhs", "rhs", "est")]
  
  # Construct edge list from coefficient table
  
  nodes <- unique(as.character(c(b$lhs, b$rhs)))
  to <- factor(b$lhs, levels=nodes) 
  from <- factor(b$rhs, levels=nodes)
  
  # Simplify node labels
  
  to <- factor(node.names[to], levels=node.names) 
  from <- factor(node.names[from], levels=node.names)
  
  return(data.frame(to, from, b=b$est))
}

nodes <- c("Phyto", "Zoo", "Fish", "Long", "TP", "TOC")
SEM1.network <- sem.network(sem1, node.names=nodes)
SEM2.network <- sem.network(sem2, node.names=nodes)

SEM1.graph <- with(SEM1.network, sem.graph(to, from, b))
SEM2.graph <- with(SEM2.network, sem.graph(to, from, b))

render_graph(SEM1.graph)
export_graph(SEM1.graph, file_name="Figure 4C.pdf")

render_graph(SEM2.graph)
export_graph(SEM2.graph, file_name="Figure 4D.pdf")

# Invalid asm.js: Function definition doesn't match use
# https://github.com/rich-iannone/DiagrammeRsvg/issues/7
# Seems to produce correct pdfs though...
