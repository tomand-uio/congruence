

# Plot 3-way variance partitioning between environment (E), space (S), 
# biotic interactions (B), and all interactions (ES, EB, SB, ESB)
# in a circular diagram with a central area representing unresolved
# variance (ESB), and outer rim representing total explained variance,
# and an annulus with color coded sectors representing variance explained by
# E, ES, S, SB, B, BE. Areas are scaled to their contributions to total variance

circle <- function(r=1) {
  v <- seq(0, 2 * pi, pi / 200)
  v <- c(v, 0)
  
  x <- r * cos(v)
  y <- r * sin(v)
  
  return(list(x=x,y=y))
}

sector <- function(r1, r2, v1, v2) {
  v <- seq(v1, v2, pi / 200)
  
  x <- c(r2 * cos(v), r1 * rev(cos(v)), r2 * cos(v1))
  y <- c(r2 * sin(v), r1 * rev(sin(v)), r2 * sin(v1))
  
  return(list(x=x,y=y))
}

# The expected order of variance partitions (V) is: E, S, B, ES, EB, SB, ESB)
# Square root-scaling radii means areas are proportional to V components

wheel.diagram <- function(V, main) {
  
  r1 <- V[7] # 111
  r2 <- sum(V)
  
  plot(c(-1,1), c(-1,1), asp=1, axes=FALSE, type="n", xlab="", ylab="",
       main=paste(main, "\n", round(100 * r2, digits=1), "%"))
  
  v <- 2 * pi * cumsum(c(0, V[c(1,4,2,6,3,5)])) / sum(V[1:6])
  
  v.names <- c("E", "ES", "S", "SB", "B", "EB")
  v.col <- c(rgb(1,0,0,0.3), rgb(1,1,0,0.3), rgb(0,1,0,0.3), 
             rgb(0,1,1,0.3), rgb(0,0,1,0.3), rgb(1,0,1,0.3))
  
  for (i in 1:6) {		
    polygon(sector(sqrt(r1), sqrt(r2), v[i], v[i+1]), 
            col=v.col[i], border=NA)
    text(sqrt((1/2) * (r1 + r2)) * cos((v[i] + v[i+1]) / 2),
         sqrt((1/2) * (r1 + r2)) * sin((v[i] + v[i+1]) / 2), 
         v.names[i], col=gray(0, 0.5))
  }
  
  polygon(circle(1))
  polygon(circle(sqrt(r2)), border=gray(0.5,0.5), lwd=3)
  polygon(circle(sqrt(r1)), col=gray(0.5,0.5), border=NA)
  text(0, 0, "ESB", col=gray(1, 0.5))
}
